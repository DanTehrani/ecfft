use crate::curve::{ec_add, ec_mul, AffinePoint, GoodCurve};
use halo2curves::ff::PrimeField;
use num_bigint::BigUint;

// 2 x 2 matrix
#[derive(Debug, Clone, Copy)]
pub struct Matrix2x2<F: PrimeField> {
    pub a: F,
    pub b: F,
    pub c: F,
    pub d: F,
}

impl<F: PrimeField> Matrix2x2<F> {
    pub fn invert(&self) -> Self {
        let det = self.a * self.d - self.b * self.c;
        assert_ne!(det, F::ZERO);
        let inv_det = det.invert().unwrap();

        Self {
            a: self.d * inv_det,
            b: -self.b * inv_det,
            c: -self.c * inv_det,
            d: self.a * inv_det,
        }
    }

    pub fn mul_vec(&self, v: (F, F)) -> (F, F) {
        let r1 = self.a * v.0 + self.b * v.1;
        let r2 = self.c * v.0 + self.d * v.1;

        (r1, r2)
    }
}

pub fn prepare_matrices<F: PrimeField>(
    L: &Vec<Vec<F>>,
) -> (Vec<Vec<Matrix2x2<F>>>, Vec<Vec<Matrix2x2<F>>>) {
    let k = L.len();
    let mut matrices = Vec::with_capacity(k);
    let mut inverse_matrices = Vec::with_capacity(k);

    for i in 0..k {
        let L_i = &L[i];
        let s_i = L_i.iter().step_by(2).map(|x| *x).collect::<Vec<F>>();
        let s_prime_i = L_i
            .iter()
            .skip(1)
            .step_by(2)
            .map(|x| *x)
            .collect::<Vec<F>>();

        let n = s_i.len();
        let nn = n / 2;
        // Matrices to go from P to P0 , P1 over the domain s
        let mut inverse_matrices_i = Vec::with_capacity(nn);
        // Matrices to go from P0, P1 to P over the extended domain s'
        let mut matrices_i = Vec::with_capacity(nn);

        for i in 0..(nn) {
            let s0 = s_i[i];
            let s1 = s_i[i + nn];

            // The denominator v(x) of the degree-2 map is x
            let q = nn - 1;

            let a = s0.pow(&[q as u64, 0, 0, 0]);
            let b = s0 * a;
            let c = s1.pow(&[q as u64, 0, 0, 0]);
            let d = s1 * c;

            let m = Matrix2x2 { a, b, c, d };
            inverse_matrices_i.push(m.invert());
        }
        inverse_matrices.push(inverse_matrices_i);

        for i in 0..(nn) {
            let s0 = s_prime_i[i];
            let s1 = s_prime_i[i + nn];

            // The denominator v(x) of the degree-2 map is x
            let q = nn - 1;

            let a = s0.pow(&[q as u64, 0, 0, 0]);
            let b = s0 * a;
            let c = s1.pow(&[q as u64, 0, 0, 0]);
            let d = s1 * c;

            let m = Matrix2x2 { a, b, c, d };
            matrices_i.push(m);
        }

        matrices.push(matrices_i);
    }

    (matrices, inverse_matrices)
}

fn to_unique<F: PrimeField>(a: &[F]) -> Vec<F> {
    let mut unique = vec![];
    for a_i in a {
        if !unique.contains(a_i) {
            unique.push(*a_i);
        }
    }
    unique
}

fn isogeny_chain<F: PrimeField>(good_curve: GoodCurve<F>) -> Vec<GoodCurve<F>> {
    let mut curves = Vec::with_capacity(good_curve.k - 1);
    curves.push(good_curve);
    for i in 0..(good_curve.k - 2) {
        curves.push(curves[i].good_iso());
    }
    curves
}

pub fn prepare_domain<F: PrimeField>(good_curve: GoodCurve<F>) -> Vec<Vec<F>> {
    let curves = isogeny_chain(good_curve.clone());

    let g = AffinePoint::new(good_curve.gx, good_curve.gy);
    // Generate the group of order 2^k and collect the x-coordinates of the points

    let G = (0..2u32.pow(good_curve.k as u32))
        .map(|i| {
            let s = BigUint::from(i);
            ec_mul::<F, GoodCurve<F>>(&g, &s, &good_curve)
        })
        .collect::<Vec<AffinePoint<F>>>();

    let mut rng = rand::thread_rng();
    let mut coset_offset_x = F::ZERO;
    let mut coset_offset_y = F::ZERO;

    loop {
        coset_offset_x = F::random(&mut rng);
        let y = (coset_offset_x
            * (coset_offset_x.square() + curves[0].a * coset_offset_x + curves[0].B_sqrt.square()))
        .sqrt();

        if y.is_some().into() {
            coset_offset_y = y.unwrap();
            break;
        }
    }

    let coset_offset = AffinePoint::new(coset_offset_x, coset_offset_y);

    // Apply the coset offset and get the projection on F
    let L0 = G
        .iter()
        .map(|p| ec_add::<F, GoodCurve<F>>(p, &coset_offset, &good_curve).x)
        .collect::<Vec<F>>();

    let mut L = Vec::with_capacity(good_curve.k);
    L.push(L0);

    // TODO: Stop doing to_unique here
    for i in 0..(good_curve.k - 2) {
        let L_i = &L[i];
        L.push(to_unique(
            &L_i.iter()
                .map(|L_i| curves[i].iso_x(*L_i))
                .collect::<Vec<F>>(),
        ))
    }

    L
}

#[cfg(test)]
mod tests {
    use super::*;
    use halo2curves::secp256k1::Fp;
    #[test]
    fn test_prepare_matrices() {
        let k = 7;
        let good_curve = GoodCurve::<Fp>::find_k(k);
        let domain = prepare_domain(good_curve);
        let (matrices, inverse_matrices) = prepare_matrices(&domain);
        for i in 0..matrices.len() {
            assert_eq!(matrices[i].len(), 2usize.pow((k - 2 - i) as u32));
        }
    }

    #[test]
    fn test_iso_chain() {
        let good_curve = GoodCurve::<Fp>::find_k(7);
        let curves = isogeny_chain(good_curve.clone());

        // The generator of the group of order 4 should be as expected
        let order_4_gx = curves[curves.len() - 1];
        let B_sqrt = order_4_gx.B_sqrt;
        assert!(order_4_gx.gx == B_sqrt || order_4_gx.gx == -B_sqrt);
    }

    #[test]
    fn test_prepare_domain() {
        let good_curve = GoodCurve::<Fp>::find_k(7);
        let domain = prepare_domain(good_curve);
        for i in 0..domain.len() {
            assert_eq!(
                domain[i].len(),
                2usize.pow((good_curve.k - i) as u32),
                "{}",
                i
            );
        }
    }
}
