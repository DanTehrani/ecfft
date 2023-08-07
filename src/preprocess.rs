use crate::curve::{ec_add, AffinePoint, GoodCurve};
use halo2curves::group::ff::PrimeField;

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
        assert_ne!(det, F::zero());
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

            let a = s0.pow_vartime(&[q as u64, 0, 0, 0]);
            let b = s0 * a;
            let c = s1.pow_vartime(&[q as u64, 0, 0, 0]);
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

            let a = s0.pow_vartime(&[q as u64, 0, 0, 0]);
            let b = s0 * a;
            let c = s1.pow_vartime(&[q as u64, 0, 0, 0]);
            let d = s1 * c;

            let m = Matrix2x2 { a, b, c, d };
            matrices_i.push(m);
        }

        matrices.push(matrices_i);
    }

    (matrices, inverse_matrices)
}

fn isogeny_chain<F: PrimeField>(good_curve: GoodCurve<F>) -> Vec<GoodCurve<F>> {
    let mut curves = Vec::with_capacity(good_curve.k - 1);
    curves.push(good_curve);
    for i in 0..(good_curve.k - 2) {
        curves.push(curves[i].good_iso());
    }
    curves
}

pub fn prepare_domain<F: PrimeField>(
    good_curve: GoodCurve<F>,
    coset_offset_x: F,
    coset_offset_y: F,
) -> Vec<Vec<F>> {
    let curves = isogeny_chain(good_curve.clone());

    let coset_offset = AffinePoint::new(coset_offset_x, coset_offset_y);
    let g = AffinePoint::new(good_curve.gx, good_curve.gy);
    let coset_g = ec_add::<F, GoodCurve<F>>(&g, &coset_offset, &good_curve);
    // Generate the group of order 2^k and collect the x-coordinates of the points

    let n = 2usize.pow(good_curve.k as u32);
    let mut G = Vec::with_capacity(n);
    G.push(coset_offset);
    G.push(coset_g.clone());

    for i in 1..n - 1 {
        G.push(ec_add::<F, GoodCurve<F>>(&G[i as usize], &g, &curves[0]))
    }

    // Collect the x-coordinates of the points
    let L0 = G.iter().map(|p| p.x).collect::<Vec<F>>();

    let mut L = Vec::with_capacity(good_curve.k);
    L.push(L0);

    for i in 0..(good_curve.k - 2) {
        let L_i = &L[i];
        L.push(
            L_i[..(L_i.len() / 2)]
                .iter()
                .map(|L_i| curves[i].iso_x(*L_i))
                .collect::<Vec<F>>(),
        )
    }

    L
}

#[cfg(test)]
mod tests {
    use crate::utils::find_coset_offset;

    use super::*;
    use halo2curves::secp256k1::Fp;
    #[test]
    fn test_prepare_matrices() {
        let k = 7;
        let good_curve = GoodCurve::<Fp>::find_k(k);
        let (coset_offset_x, coset_offset_y) =
            find_coset_offset(good_curve.a, good_curve.B_sqrt.square());
        let domain = prepare_domain(good_curve, coset_offset_x, coset_offset_y);
        let (matrices, _) = prepare_matrices(&domain);
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
        let good_curve = GoodCurve::<Fp>::find_k(5);
        let (coset_offset_x, coset_offset_y) =
            find_coset_offset(good_curve.a, good_curve.B_sqrt.square());
        let domain = prepare_domain(good_curve.clone(), coset_offset_x, coset_offset_y);
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
