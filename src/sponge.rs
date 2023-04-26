use crate::poseidon::k256_consts;
use crate::poseidon::{Poseidon, PoseidonConstants};
use ff::PrimeField;
use sha3::{Digest, Sha3_256};
use std::result::Result;

#[derive(Clone)]
pub enum SpongeOp {
    Absorb(usize),
    Squeeze(usize),
}

#[derive(Clone)]
pub struct IOPattern(pub Vec<SpongeOp>);

// Implements SAFE (Sponge API for Field Elements): https://hackmd.io/bHgsH6mMStCVibM_wYvb2w
pub struct PoseidonSponge<F: PrimeField> {
    pub absorb_pos: usize,
    pub squeeze_pos: usize,
    pub io_count: usize,
    pub io_pattern: Option<IOPattern>,
    pub rate: usize,
    pub capacity: usize,
    poseidon: Poseidon<F>,
}

pub enum SpongeCurve {
    K256,
}

impl<F: PrimeField<Repr = [u8; 32]>> PoseidonSponge<F> {
    pub fn construct(
        domain_separator: &[u8],
        curve: SpongeCurve,
        io_pattern: Option<IOPattern>,
    ) -> Self {
        // Parse the constants from string
        let constants = match curve {
            SpongeCurve::K256 => {
                let round_constants: Vec<F> = k256_consts::ROUND_CONSTANTS
                    .iter()
                    .map(|x| F::from_str_vartime(x).unwrap())
                    .collect();

                let mds_matrix: Vec<Vec<F>> = k256_consts::MDS_MATRIX
                    .iter()
                    .map(|x| {
                        x.iter()
                            .map(|y| F::from_str_vartime(y).unwrap())
                            .collect::<Vec<F>>()
                    })
                    .collect();

                PoseidonConstants::new(
                    round_constants,
                    mds_matrix,
                    k256_consts::NUM_FULL_ROUNDS,
                    k256_consts::NUM_PARTIAL_ROUNDS,
                )
            }
        };

        let tag = Self::compute_tag(domain_separator, &io_pattern);

        let state = vec![tag, F::zero(), F::zero()];

        let poseidon = Poseidon::new(constants, state);

        Self {
            absorb_pos: 0,
            squeeze_pos: 0,
            io_count: 0,
            io_pattern,
            rate: 2,
            capacity: 1,
            poseidon,
        }
    }

    // Compute tag as described in section 2.3 of the SAFE documentation
    fn compute_tag(domain_separator: &[u8], io_pattern: &Option<IOPattern>) -> F {
        // step 1: Encode
        let io_words = match &io_pattern {
            Some(io) => {
                io.0.iter()
                    .map(|io_i| match io_i {
                        SpongeOp::Absorb(n) => (n + 0x80000000) as u32,
                        SpongeOp::Squeeze(n) => (*n) as u32,
                    })
                    .collect()
            }
            None => {
                vec![]
            }
        };

        // step 2: Aggregate
        let mut io_words_aggregated = vec![];
        for io_word in io_words {
            if io_words_aggregated.len() == 0 {
                io_words_aggregated.push(io_word);
            } else {
                let i = io_words_aggregated.len() - 1;
                if io_words_aggregated[i] > 0x80000000 && io_word > 0x80000000 {
                    io_words_aggregated[i] += io_word - 0x80000000;
                } else if io_words_aggregated[i] < 0x80000000 && io_word < 0x80000000 {
                    io_words_aggregated[i] += io_word;
                } else {
                    io_words_aggregated.push(io_word);
                }
            }
        }

        // step 3: Serialize
        let mut io_bytes = vec![];
        for io_word in io_words_aggregated {
            io_word.to_be_bytes().iter().for_each(|x| io_bytes.push(*x));
        }
        io_bytes.extend_from_slice(domain_separator);

        // step 4: Hash
        let mut hasher = Sha3_256::new();
        hasher.update(&io_bytes.as_slice());
        let result = hasher.finalize();

        // Truncate the first 128 bits of the hash to compute the tag
        let mut tag = Vec::with_capacity(32);
        tag.extend_from_slice(&[0u8; 16]);
        tag.extend_from_slice(&result[0..16]);

        F::from_repr(tag.as_slice().try_into().unwrap()).unwrap()
    }

    pub fn absorb(&mut self, x: &[F]) {
        if x.len() == 0 {
            return;
        }

        for x_i in x {
            if self.absorb_pos == self.rate {
                self.permute();
                self.absorb_pos = 0
            }

            self.poseidon.state[self.absorb_pos] = *x_i;
            self.absorb_pos += 1;
        }

        // TODO: Verify the IO pattern
        self.io_count += 1;
        self.squeeze_pos = self.rate;
    }

    pub fn squeeze(&mut self, length: usize) -> Vec<F> {
        let mut y = Vec::with_capacity(length);
        if length == 0 {
            return vec![];
        }

        for _ in 0..length {
            if self.squeeze_pos == self.rate {
                self.permute();
                self.squeeze_pos = 0;
                self.absorb_pos = 0;
            }

            y.push(self.poseidon.state[self.squeeze_pos]);
            self.squeeze_pos += 1;
        }

        self.io_count += 1;
        y
    }

    pub fn finish(&self) -> Result<(), String> {
        match self.io_pattern {
            None => return Ok(()),
            Some(ref io_pattern) => {
                if self.io_count != io_pattern.0.len() {
                    return Err("IO pattern mismatch".to_string());
                }
            }
        }

        Ok(())
    }

    fn permute(&mut self) {
        self.poseidon.permute();
        self.poseidon.pos = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //    use secq256k1::field::field_secq::FieldElement as Fp;
    use halo2curves::secp256k1::Fp;

    #[test]
    fn test_interactive_protocol() {
        let io_pattern = IOPattern(vec![
            SpongeOp::Absorb(2),
            SpongeOp::Squeeze(1),
            SpongeOp::Absorb(1),
            SpongeOp::Squeeze(3),
        ]);

        let io = vec![vec![Fp::from(1), Fp::from(2)], vec![Fp::from(3)]].concat();

        let mut sponge =
            PoseidonSponge::construct(b"test", SpongeCurve::K256, Some(io_pattern.clone()));

        let mut io_position = 0;
        for op in io_pattern.0 {
            match op {
                SpongeOp::Absorb(l) => {
                    sponge.absorb(&io[io_position..(io_position + l)]);
                    io_position += l;
                }
                SpongeOp::Squeeze(l) => {
                    sponge.squeeze(l);
                }
            }
        }

        assert_eq!(sponge.finish(), Ok(()));
    }
}
