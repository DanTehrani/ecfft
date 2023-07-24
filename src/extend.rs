use crate::curve::GoodCurve;
use ff::PrimeField;

// https://solvable.group/posts/ecfft/

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

        if det.is_zero().into() {
            return self.clone();
        }

        let inv_det = det.invert().unwrap();
        debug_assert!(det != F::ZERO);

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

// L: Evaluation domain
pub fn compute_matrices<F: PrimeField>(
    L: &[F],
    curve: &GoodCurve<F>,
) -> (Vec<Vec<Matrix2x2<F>>>, Vec<Vec<Matrix2x2<F>>>) {
    let mut n = L.len();

    let log_n = (L.len() as f64).log2() as usize;
    let mut matrices = Vec::with_capacity(log_n);
    let mut inverse_matrices = Vec::with_capacity(log_n);

    for i in 0..log_n {
        let nn = n / 2;
        let mut matrices_i = Vec::with_capacity(nn);
        let mut inverse_matrices_i = Vec::with_capacity(nn);

        for i in 0..(nn) {
            let s0 = L[i];
            let s1 = L[i + nn];

            // The denominator v(x) of the degree-2 map is x
            let q = nn - 1;

            let a = s0.pow(&[q as u64, 0, 0, 0]);
            let b = s0 * a;
            let c = s1.pow(&[q as u64, 0, 0, 0]);
            let d = s1 * c;

            let m = Matrix2x2 { a, b, c, d };
            matrices_i.push(m);
            inverse_matrices_i.push(m.invert());
        }

        matrices.push(matrices_i);
        inverse_matrices.push(inverse_matrices_i);
        n >>= 1
    }

    (matrices, inverse_matrices)
}

pub fn compute_domain<F: PrimeField>() -> Vec<Vec<F>> {
    todo!();
}

// Low-degree extension
pub fn extend<F: PrimeField>(
    evals: &[F],
    L: &Vec<Vec<F>>,
    matrices: &Vec<Vec<Matrix2x2<F>>>,
    inverse_matrices: &Vec<Vec<Matrix2x2<F>>>,
    i: usize,
) -> Vec<F> {
    if evals.len() == 1 {
        return evals.to_vec();
    }

    let s = L[i].iter().step_by(2).map(|x| *x).collect::<Vec<F>>();
    let s_prime = L[i]
        .iter()
        .skip(1)
        .step_by(2)
        .map(|x| *x)
        .collect::<Vec<F>>();

    debug_assert_eq!(s.len(), s_prime.len());
    debug_assert_eq!(evals.len(), s.len());
    let n = s.len();
    let nn = n / 2;

    if s.len() == 1 {
        return evals.to_vec();
    }

    // Deduce the evaluation of Q_1, Q_2, over \psi(s0)
    // from the evaluation Q over s0 and s1
    let mut Q_1_evals = vec![];
    let mut Q_2_evals = vec![];

    for j in 0..nn {
        let y0 = evals[j];
        let y1 = evals[j + nn];

        let m_i = matrices[i][j];
        // Evaluations of Q_0, Q_1 at t
        let (q0, q1) = m_i.mul_vec((y0, y1));

        Q_1_evals.push(q0);
        Q_2_evals.push(q1);
    }

    let Q_1_evals_prime = extend(&Q_1_evals, L, matrices, inverse_matrices, i + 1);
    let Q_2_evals_prime = extend(&Q_2_evals, L, matrices, inverse_matrices, i + 1);

    let mut extended_evals = vec![F::ZERO; evals.len() * 2];

    for (i, ((q1, q2), m)) in Q_1_evals_prime
        .iter()
        .zip(Q_2_evals_prime.iter())
        .zip(inverse_matrices[i].iter())
        .enumerate()
    {
        let (e1, e2) = m.mul_vec((*q1, *q2));
        extended_evals[i] = e1;
        extended_evals[i + nn] = e2;
    }

    extended_evals
}

#[cfg(test)]
mod tests {
    use crate::curve::{ec_add, ec_mul};
    use crate::Fp;
    use ff::{Field, PrimeField};
    use num_bigint::BigUint;

    use super::*;

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
    fn test_extend() {
        // 0x0afd2c697fc6654762108112f46adc57a7d7ccf815da9e8dc41612a308e5a692
        let a = Fp::from_str_vartime(
            "14757638988853573751207811144824545048510466205597221722951176097150491841035",
        )
        .unwrap();

        let sqrt_b = Fp::from_str_vartime(
            "5567257267101488941029590412198289470015986673315929203017157096142767142896",
        )
        .unwrap();

        let gx = Fp::from_str_vartime(
            "21255228318590632286479157202391852408063197350582547415436433972957685992507",
        )
        .unwrap();
        let gy = Fp::from_str_vartime(
            "26972981084054348264666231698996206396459155228721448282357884867841822592398",
        )
        .unwrap();

        let g = Point { x: gx, y: gy };

        let curve_k = 4;
        // let curve = GoodCurve::new(a, sqrt_b, g, curve_k);
        let curves = good_isogeny_chain(curve_k);
        let curve = curves[0];

        let G = (0..2u32.pow(curve.k as u32))
            .map(|i| {
                let s = BigUint::from(i);
                ec_mul::<F, GoodCurve<F>>(&curve.two_sylow_generator, &s, &curve)
            })
            .collect::<Vec<Point<F>>>();

        let coset_offset =
           // ec_mul::<F, GoodCurve<F>>(&curve.two_sylow_generator, &BigUint::from(100u32), &curve);
           Point::new(F::from_str_vartime("105623886150579165427389078198493427091405550492761682382732004625374789850161").unwrap(), F::from_str_vartime("7709812624542158994629670452026922591039826164720902911013234773380889499231").unwrap());

        let log_n = curve.k;
        let n = 2usize.pow(log_n as u32);

        let L = G
            .iter()
            //  .map(|g| ec_add::<F, GoodCurve<F>>(&g, &coset_offset, &curve).x)
            .map(|g| g.x)
            .collect::<Vec<F>>();

        let (matrices, inverse_matrices) = compute_matrices(&L, &curve);

        let mut L_all = Vec::with_capacity(log_n);
        L_all.push(L.clone());

        for i in 0..(log_n - 1) {
            let L_i = &L_all[i];
            L_all.push(to_unique(
                &L_i.iter()
                    .map(|L_i| curves[i].iso_x(*L_i))
                    .collect::<Vec<F>>(),
            ))
        }

        let coeffs = (0..(L.len() / 2))
            .map(|i: usize| F::from(i as u64))
            .collect::<Vec<F>>();

        let poly = UniPoly::new(coeffs);

        let s = L.iter().step_by(2).map(|x| *x);
        // Evaluations of poly over s
        let evals = s.map(|s_i| poly.eval(s_i)).collect::<Vec<F>>();

        // Evaluations of poly over L
        let extended_evals = encode(&evals, &L_all, &matrices, &inverse_matrices);

        println!("extended_evals {:?}", extended_evals);
        // let expected_evals = L.iter().map(|L_i| poly.eval(*L_i)).collect::<Vec<F>>();
        // assert_eq!(extended_evals, expected_evals);
    }
}
