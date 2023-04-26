use crate::poseidon::k256_consts;
use crate::poseidon::k256_consts::{NUM_FULL_ROUNDS, NUM_PARTIAL_ROUNDS};
use crate::poseidon::{Poseidon, PoseidonConstants};
use ff::PrimeField;
use std::result::Result;

#[derive(Clone)]
pub enum SpongeOp {
    Absorb(usize),
    Squeeze(usize),
}

#[derive(Clone)]
pub struct IOPattern(pub Vec<SpongeOp>);

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

impl<F: PrimeField> PoseidonSponge<F> {
    pub fn construct(curve: SpongeCurve, io_pattern: Option<IOPattern>) -> Self {
        // Compute the tag T
        // Set the permutation state to all zeros and add T
        // to the first 128 bits of the state

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
        let poseidon = Poseidon::new(constants);

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

    pub fn absorb(&mut self, x: &[F]) {
        if x.len() == 0 {
            return;
        }

        for (i, x_i) in x.iter().enumerate() {
            if self.absorb_pos == self.rate {
                self.permute();
                self.absorb_pos = 0
            }

            self.poseidon.state[self.absorb_pos] = *x_i;
            self.absorb_pos += 1;
        }

        // Verify the IO pattern
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

    fn permute(&mut self) {
        self.poseidon.permute();
        self.poseidon.pos = 0;
    }

    fn finish(&self) -> Result<(), String> {
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

        let mut sponge = PoseidonSponge::construct(SpongeCurve::K256, Some(io_pattern.clone()));

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
