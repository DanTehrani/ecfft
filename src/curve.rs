use ark_ff::PrimeField;

// Copied from https://github.com/andrewmilson/ecfft/blob/main/src/ec.rs
use std::fmt::Debug;

use crate::{get_2sylow_subgroup::det_2sylow_cyclic_subgroup, utils::is_quad_residue};

#[derive(Debug, Clone, Copy)]
pub struct AffinePoint<F: PrimeField> {
    pub x: F,
    pub y: F,
    pub infinity: bool,
}

impl<F: PrimeField> AffinePoint<F> {
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }

    // Point at infinity
    pub fn infinity() -> Self {
        Self {
            x: F::zero(),
            y: F::zero(),
            infinity: true,
        }
    }
}

// Good curve (Definition 4.1 https://www.math.toronto.edu/swastik/ECFFT2.pdf)
// over an odd size field
#[derive(Clone, Copy, Debug)]
pub struct GoodCurve<F: PrimeField> {
    pub a: F,
    pub B_sqrt: F, // Sqrt of the B coefficient
    pub gx: F,     // Generator of the 2-Sylow subgroup
    pub gy: F,     // Generator of the 2-Sylow subgroup
    pub k: usize,  // 2^k = n is the order of the curve
}

impl<F: PrimeField> GoodCurve<F> {
    pub fn new(a: F, B_sqrt: F, gx: F, gy: F, k: usize) -> Self {
        debug_assert!(B_sqrt != F::zero());
        let b = B_sqrt.square();

        // Check a - 2B != 0 where B is the degree 1 coefficient of the curve
        debug_assert!((a.square() - b.double().double()) != F::zero());

        // Check a + 2b is a quadratic residue
        debug_assert!(is_quad_residue(a + B_sqrt.double()));

        Self {
            a,
            B_sqrt,
            gx,
            gy,
            k,
        }
    }

    pub fn find_k(k: usize) -> Self {
        // Find a curve with a 2-Sylow subgroup of order 10
        let mut rng = rand::thread_rng();
        loop {
            let a = F::rand(&mut rng);
            let sqrt_b = F::rand(&mut rng);
            match det_2sylow_cyclic_subgroup(a, sqrt_b) {
                Some((c_k, gx)) => {
                    if c_k == k {
                        let b = sqrt_b.square();
                        let gy = (gx * (gx.square() + a * gx + b)).sqrt().unwrap();
                        return Self::new(a, sqrt_b, gx, gy, k);
                    }
                }
                None => continue,
            }
        }
    }

    pub fn good_point(&self) -> AffinePoint<F> {
        let y = self.B_sqrt * (self.a + self.B_sqrt.double()).sqrt().unwrap();
        AffinePoint::new(self.B_sqrt, y)
    }

    pub fn iso_x(&self, x: F) -> F {
        if x.is_zero().into() {
            return F::zero();
        }

        (x - self.B_sqrt).square() * x.inverse().unwrap()
    }

    // Return the curve that corresponds to the codomain of a 2-isogeny
    pub fn good_iso(&self) -> Self {
        let a = self.a;
        let b_sqrt = self.B_sqrt;
        let b = b_sqrt.square();
        let g = AffinePoint::new(self.gx, self.gy);

        let a_prime = a + b_sqrt.double().double() + b_sqrt.double();
        let b_prime = (a * b_sqrt).double().double() + F::from(8u64) * b;

        let b_sqrt_prime = b_prime.sqrt();
        if b_sqrt_prime.is_none().into() {
            panic!(
                "b_prime is not a square b_prime {:?} {:?}",
                b_prime, b_sqrt_prime
            );
        }

        let gx_prime = self.iso_x(g.x);
        let gy_prime = (gx_prime * (gx_prime.square() + a_prime * gx_prime + b_prime))
            .sqrt()
            .unwrap();
        let k_prime = self.k - 1;

        let E_prime = GoodCurve::new(a_prime, b_sqrt_prime.unwrap(), gx_prime, gy_prime, k_prime);
        E_prime
    }
}

/// General Weierstrass curve of the form:
/// `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`
pub trait WeierstrassCurve: Clone + Copy + Debug {
    type F: PrimeField;

    /// Returns Weierstrass equation coefficient `a1`
    fn a1(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a2`
    fn a2(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a3`
    fn a3(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a4`
    fn a4(&self) -> Self::F;

    /// Returns Weierstrass equation coefficient `a6`
    fn a6(&self) -> Self::F;
}

impl<F: PrimeField> WeierstrassCurve for GoodCurve<F> {
    type F = F;

    fn a1(&self) -> Self::F {
        F::zero()
    }

    fn a2(&self) -> Self::F {
        self.a
    }

    fn a3(&self) -> Self::F {
        F::zero()
    }

    fn a4(&self) -> Self::F {
        self.B_sqrt.square()
    }

    fn a6(&self) -> Self::F {
        F::zero()
    }
}

pub fn ec_add<F: PrimeField, C: WeierstrassCurve>(
    a: &AffinePoint<F>,
    b: &AffinePoint<F>,
    curve: &GoodCurve<F>,
) -> AffinePoint<F> {
    // Addition law for elliptic curve groups
    // Source: "The arithmetic of elliptic curves, 2nd ed." Silverman, III.2.3
    if a.infinity {
        *b
    } else if b.infinity {
        *a
    } else {
        let a1 = curve.a1();
        let a2 = curve.a2();
        let a3 = curve.a3();
        let a4 = curve.a4();
        let a6 = curve.a6();
        let x1 = a.x;
        let y1 = a.y;
        let x2 = b.x;
        let y2 = b.y;

        if x1 == x2 && (y1 + y2 + a1 * x2 + a3).is_zero().into() {
            AffinePoint::infinity()
        } else {
            let lambda: F;
            let nu: F;
            if x1 == x2 {
                // tangent line
                let x1x1 = x1.square();
                let a2x1 = a2 * x1;
                let a1x1 = a1 * x1;
                lambda = (x1x1 + x1x1 + x1x1 + a2x1 + a2x1 + a4 - a1 * y1)
                    * (y1 + y1 + a1x1 + a3).inverse().unwrap();
                nu = (-(x1x1 * x1) + a4 * x1 + a6 + a6 - a3 * y1)
                    * (y1 + y1 + a1 * x1 + a3).inverse().unwrap();
            } else {
                // slope through the AffinePoints
                lambda = (y2 - y1) * (x2 - x1).inverse().unwrap();
                nu = (y1 * x2 - y2 * x1) * (x2 - x1).inverse().unwrap();
            }
            let x3 = lambda.square() + a1 * lambda - a2 - x1 - x2;
            let y3 = -(lambda + a1) * x3 - nu - a3;
            AffinePoint::new(x3, y3)
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::find_coset_offset;

    use super::*;
    use ark_secp256k1::Fq as Fp;

    #[test]
    fn test_find_k() {
        let k = 9;
        let curve = GoodCurve::<Fp>::find_k(k);

        assert_eq!(curve.k, k);

        // Should generate a group of order 2^k
        let g = AffinePoint::new(curve.gx, curve.gy);
        let n = 2usize.pow(k as u32);
        let mut G = Vec::with_capacity(n);
        G.push(AffinePoint::infinity());
        G.push(g.clone());

        for i in 1..2usize.pow(k as u32) - 1 {
            G.push(ec_add::<Fp, GoodCurve<Fp>>(&G[i], &g, &curve));
        }

        let (coset_offset_x, coset_offset_y) =
            find_coset_offset(curve.a, curve.B_sqrt * curve.B_sqrt);
        let coset_offset = AffinePoint::new(coset_offset_x, coset_offset_y);

        let L = G
            .iter()
            .map(|G_i| ec_add::<Fp, GoodCurve<Fp>>(G_i, &coset_offset, &curve).x)
            .collect::<Vec<Fp>>();

        assert_eq!(L.len(), 2usize.pow(k as u32));
    }
}
