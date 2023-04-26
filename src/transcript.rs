use crate::sponge::{PoseidonSponge, SpongeCurve};
use ff::PrimeField;
use halo2curves::{CurveAffineExt, FieldExt};

pub struct PoseidonTranscript<C: CurveAffineExt> {
    sponge: PoseidonSponge<C::ScalarExt>,
}

impl<C> PoseidonTranscript<C>
where
    C: CurveAffineExt,
    C::ScalarExt: FieldExt<Repr = [u8; 32]>,
    C::Base: FieldExt<Repr = [u8; 32]>,
{
    pub fn new(domain_separator: &[u8], curve: SpongeCurve) -> Self {
        // The scalar field of the curve specified by the generic argument
        // is used as the finite field of the Poseidon sponge.
        match curve {
            SpongeCurve::K256 => {
                assert_eq!(
                    C::ScalarExt::MODULUS,
                    "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f"
                );
            }
        }

        Self {
            sponge: PoseidonSponge::construct(domain_separator, curve, None),
        }
    }

    // Append a byte array to the transcript.
    pub fn append_bytes(&mut self, bytes: &[u8]) {
        assert!(bytes.len() <= 64);
        let mut padded_bytes = Vec::with_capacity(64);
        padded_bytes.extend_from_slice(bytes);
        padded_bytes.resize(64, 0);
        self.sponge.absorb(&[C::ScalarExt::from_bytes_wide(
            padded_bytes.as_slice().try_into().unwrap(),
        )]);
    }

    // Append a group element to the transcript.
    pub fn append_point(&mut self, point: &C) {
        let coords = point.coordinates().unwrap();
        let x = coords.x();
        let y = coords.y();

        let x: [u8; 32] = x.to_repr();
        let y = y.to_repr();

        self.append_bytes(&x);
        self.append_bytes(&y);
    }

    // Append a scalar field element to the transcript.
    pub fn append_scalar(&mut self, fe: &C::ScalarExt) {
        self.sponge.absorb(&[*fe]);
    }

    // Squeeze a vector of scalar field elements from the transcript.
    pub fn squeeze(&mut self, length: usize) -> Vec<C::ScalarExt> {
        self.sponge.squeeze(length)
    }
}
