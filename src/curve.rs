use ff::PrimeField;
use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::Zero;

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
            x: F::ZERO,
            y: F::ZERO,
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
        assert!(B_sqrt != F::ZERO);
        let b = B_sqrt.square();

        // Check a - 2B != 0 where B is the degree 1 coefficient of the curve
        assert!((a.square() - b.double().double()) != F::ZERO);

        // Check a + 2b is a quadratic residue
        assert!(is_quad_residue(a + B_sqrt.double()));

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
            let a = F::random(&mut rng);
            let sqrt_b = F::random(&mut rng);
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

    /*
    pub fn iso(&self, p: &AffinePoint<F>) -> AffinePoint<F> {
        let b = self.B_sqrt.square();

        if p.is_zero() {
            return AffinePoint::zero();
        }

        let good_point = self.good_point();

        if good_point.x == p.x {
            return AffinePoint::zero();
        }

        // If p.x is the identity, then the isogeny is the zero map (?)
        if self.gx == p.x {
            return AffinePoint::zero();
        }

        let x_prime = self.iso_x(p.x);
        let y_prime = (F::ONE - b * p.x.square().invert().unwrap()) * p.y;

        AffinePoint {
            x: x_prime,
            y: y_prime,
        }
    }
    */

    pub fn iso_x(&self, x: F) -> F {
        if x.is_zero().into() {
            return F::ZERO;
        }

        /*
        if self.gx == x {
            return F::ZERO;
        }
        */

        (x - self.B_sqrt).square() * x.invert().unwrap()

        /*
        let b = self.B_sqrt.square();
        if x.is_zero().into() {
            return F::ZERO;
        }

        let good_point = self.good_point();

        if good_point.x == x {
            return F::ZERO;
        }

        if self.gx == x {
            return F::ZERO;
        }

        x - self.B_sqrt.double() + b * x.invert().unwrap()
        */
    }

    // Return the curve that corresponds to the codomain of a 2-isogeny
    pub fn good_iso(&self) -> Self {
        let a = self.a;
        let b_sqrt = self.B_sqrt;
        let b = b_sqrt.square();
        let g = AffinePoint::new(self.gx, self.gy);

        let a_prime = a + b_sqrt.double().double() + b_sqrt.double();
        let b_prime = (a * b_sqrt).double().double() + F::from(8) * b;

        let b_sqrt_prime = b_prime.sqrt();
        if b_sqrt_prime.is_none().into() {
            panic!("b_prime is not a square");
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
        F::ZERO
    }

    fn a2(&self) -> Self::F {
        self.a
    }

    fn a3(&self) -> Self::F {
        F::ZERO
    }

    fn a4(&self) -> Self::F {
        self.B_sqrt.square()
    }

    fn a6(&self) -> Self::F {
        F::ZERO
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
                    * (y1 + y1 + a1x1 + a3).invert().unwrap();
                nu = (-(x1x1 * x1) + a4 * x1 + a6 + a6 - a3 * y1)
                    * (y1 + y1 + a1 * x1 + a3).invert().unwrap();
            } else {
                // slope through the AffinePoints
                lambda = (y2 - y1) * (x2 - x1).invert().unwrap();
                nu = (y1 * x2 - y2 * x1) * (x2 - x1).invert().unwrap();
            }
            let x3 = lambda.square() + a1 * lambda - a2 - x1 - x2;
            let y3 = -(lambda + a1) * x3 - nu - a3;
            AffinePoint::new(x3, y3)
        }
    }
}

pub fn ec_mul<F: PrimeField, C: WeierstrassCurve>(
    a: &AffinePoint<F>,
    s: &BigUint,
    curve: &GoodCurve<F>,
) -> AffinePoint<F> {
    let mut res = AffinePoint::infinity();
    let mut acc = a.clone();
    let mut s = s.clone();
    while s.is_zero() == false {
        if s.is_odd() {
            res = ec_add::<F, C>(&acc, &res, curve);
        }
        acc = ec_add::<F, C>(&acc, &acc, curve);
        s >>= 1;
    }
    res
}

#[cfg(test)]
mod tests {
    use ff::Field;

    use super::*;
    use crate::utils::tests::Fp;

    fn to_unique<F: PrimeField>(a: &[F]) -> Vec<F> {
        let mut unique = vec![];
        for a_i in a {
            if !unique.contains(a_i) {
                unique.push(*a_i);
            }
        }
        unique
    }

    #[test]
    fn test_find_k() {
        let k = 5;
        let curve = GoodCurve::<Fp>::find_k(k);

        assert_eq!(curve.k, k);

        // Should generate a group of order 2^k
        let G0 = AffinePoint::new(curve.gx, curve.gy);
        let mut G = vec![];
        for i in 0..2usize.pow(k as u32) {
            let s = BigUint::from(i as u32);
            G.push(ec_mul::<Fp, GoodCurve<Fp>>(&G0, &s, &curve));
        }

        let mut rng = rand::thread_rng();
        let mut coset_offset_x = Fp::ZERO;
        let mut coset_offset_y = Fp::ZERO;
        loop {
            coset_offset_x = Fp::random(&mut rng);
            let y = (coset_offset_x
                * (coset_offset_x.square() + curve.a * coset_offset_x + curve.B_sqrt.square()))
            .sqrt();

            if y.is_some().into() {
                coset_offset_y = y.unwrap();
                break;
            }
        }

        let coset_offset = AffinePoint::new(coset_offset_x, coset_offset_y);

        let L = G
            .iter()
            .map(|G_i| ec_add::<Fp, GoodCurve<Fp>>(G_i, &coset_offset, &curve).x)
            .collect::<Vec<Fp>>();

        assert_eq!(L.len(), 2usize.pow(k as u32));
    }
}
