pub(crate) mod k256_consts;
use ff::PrimeField;

pub struct PoseidonConstants<F: PrimeField> {
    pub round_keys: Vec<F>,
    pub mds_matrix: Vec<Vec<F>>,
    pub num_full_rounds: usize,
    pub num_partial_rounds: usize,
}

impl<F: PrimeField> PoseidonConstants<F> {
    pub fn new(
        round_constants: Vec<F>,
        mds_matrix: Vec<Vec<F>>,
        num_full_rounds: usize,
        num_partial_rounds: usize,
    ) -> Self {
        Self {
            num_full_rounds,
            num_partial_rounds,
            mds_matrix,
            round_keys: round_constants,
        }
    }
}

pub struct Poseidon<F: PrimeField> {
    pub state: Vec<F>,
    pub constants: PoseidonConstants<F>,
    pub pos: usize,
}

impl<F: PrimeField> Poseidon<F> {
    pub fn new(constants: PoseidonConstants<F>, state: Vec<F>) -> Self {
        Self {
            state,
            constants,
            pos: 0,
        }
    }

    pub fn permute(&mut self) {
        let full_rounds_half = self.constants.num_full_rounds / 2;

        // First half of full rounds
        for _ in 0..full_rounds_half {
            self.full_round();
        }

        // Partial rounds
        for _ in 0..self.constants.num_partial_rounds {
            self.partial_round();
        }

        // Second half of full rounds
        for _ in 0..full_rounds_half {
            self.full_round();
        }
    }

    pub fn hash(&mut self, input: Vec<F>) -> F {
        // add padding
        let mut input = input.clone();

        let domain_tag = 3; // 2^arity - 1
        input.insert(0, F::from(domain_tag));

        self.state = input;
        self.permute();

        self.state[1]
    }

    fn add_constants(&mut self) {
        // Add round constants
        for i in 0..self.state.len() {
            self.state[i] += self.constants.round_keys[i + self.pos];
        }
    }

    // MDS matrix multiplication
    fn matrix_mul(&mut self) {
        let mut result = Vec::new();

        for val in self.constants.mds_matrix.iter() {
            let mut tmp = F::zero();
            for (j, element) in self.state.iter().enumerate() {
                tmp += val[j] * element
            }
            result.push(tmp)
        }

        self.state = result;
    }

    fn full_round(&mut self) {
        let t = self.state.len();
        self.add_constants();

        // S-boxes
        for i in 0..t {
            self.state[i] = self.state[i].pow_vartime(&[5, 0, 0, 0]);
        }

        self.matrix_mul();

        // Update the position of the round constants that are added
        self.pos += self.state.len();
    }

    fn partial_round(&mut self) {
        self.add_constants();

        // S-box
        self.state[0] = self.state[0].pow_vartime(&[5, 0, 0, 0]);

        self.matrix_mul();

        // Update the position of the round constants that are added
        self.pos += self.state.len();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::secp256k1::Fp;

    #[test]
    fn test_k256() {
        let input = vec![
            Fp::from_str_vartime("1234567").unwrap(),
            Fp::from_str_vartime("109987").unwrap(),
        ];

        let round_constants: Vec<Fp> = k256_consts::ROUND_CONSTANTS
            .iter()
            .map(|x| Fp::from_str_vartime(x).unwrap())
            .collect();

        let mds_matrix: Vec<Vec<Fp>> = k256_consts::MDS_MATRIX
            .iter()
            .map(|x| {
                x.iter()
                    .map(|y| Fp::from_str_vartime(y).unwrap())
                    .collect::<Vec<Fp>>()
            })
            .collect();

        let constants = PoseidonConstants::<Fp>::new(
            round_constants,
            mds_matrix,
            k256_consts::NUM_FULL_ROUNDS,
            k256_consts::NUM_PARTIAL_ROUNDS,
        );

        let state = vec![Fp::zero(); 3];
        let mut poseidon = Poseidon::new(constants, state);

        let digest = poseidon.hash(input);

        assert_eq!(
            digest,
            Fp::from_bytes(&[
                68, 120, 17, 40, 199, 247, 48, 80, 236, 89, 92, 44, 207, 217, 83, 62, 184, 194,
                173, 48, 66, 119, 238, 98, 175, 232, 78, 234, 75, 101, 229, 148
            ])
            .unwrap()
        );
    }
}
