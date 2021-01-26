use crate::pow::Pow;
use curve25519_dalek::scalar::Scalar;
use rand::{CryptoRng, RngCore};

#[derive(Debug, Clone)]
pub struct Polynom {
    pub coeffs: Vec<Scalar>,
}

impl Polynom {
    pub fn random<R: RngCore + CryptoRng>(rng: &mut R, zero_coeff: &Scalar, order: usize) -> Self {
        let mut p = Polynom {
            coeffs: Vec::with_capacity(order),
        };

        p.coeffs.push(zero_coeff.clone());

        for _ in 1..(order + 1) {
            p.coeffs.push(Scalar::random(rng));
        }

        p
    }

    pub fn at(&self, x: &Scalar) -> Scalar {
        let mut res = self.coeffs[0].clone();
        for (i, c) in self.coeffs.iter().enumerate().skip(1) {
            res += c * x.pow(i as u64);
        }
        res
    }
}

impl Default for Polynom {
    fn default() -> Self {
        Polynom {
            coeffs: Vec::with_capacity(0),
        }
    }
}
