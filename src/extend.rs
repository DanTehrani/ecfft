use crate::preprocess::Matrix2x2;
use halo2curves::ff::PrimeField;

// https://solvable.group/posts/ecfft/

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

    let L_i = &L[i];

    debug_assert_eq!(evals.len(), L_i.len() / 2);
    let n = evals.len();
    let nn = n / 2;
    debug_assert_eq!(inverse_matrices[i].len(), nn);

    if L_i.len() == 2 {
        return evals.to_vec();
    }

    // Deduce the evaluation of Q_1, Q_2, over \psi(s0)
    // from the evaluation Q over s0 and s1
    let mut Q_1_evals = vec![];
    let mut Q_2_evals = vec![];

    for j in 0..nn {
        let y0 = evals[j];
        let y1 = evals[j + nn];

        let m_i = inverse_matrices[i][j];
        // Evaluations of Q_0, Q_1 at t
        let (q0, q1) = m_i.mul_vec((y0, y1));

        Q_1_evals.push(q0);
        Q_2_evals.push(q1);
    }

    let Q_1_evals_prime = extend(&Q_1_evals, L, matrices, inverse_matrices, i + 1);
    let Q_2_evals_prime = extend(&Q_2_evals, L, matrices, inverse_matrices, i + 1);

    let mut extended_evals = vec![F::ZERO; evals.len()];

    for (i, ((q1, q2), m)) in Q_1_evals_prime
        .iter()
        .zip(Q_2_evals_prime.iter())
        .zip(matrices[i].iter())
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
    use crate::curve::GoodCurve;
    use crate::preprocess::{prepare_domain, prepare_matrices};
    use halo2curves::bn256::Fq as Fp;

    use super::*;

    #[derive(Clone)]
    pub struct UniPoly<F>
    where
        F: PrimeField,
    {
        pub coeffs: Vec<F>,
    }

    impl<F> UniPoly<F>
    where
        F: PrimeField,
    {
        pub fn new(coeffs: Vec<F>) -> Self {
            Self { coeffs } // [x^0, x^1, x^2, x^3...]
        }

        pub fn eval(&self, x: F) -> F {
            let mut result = F::ZERO;
            for (i, coeff) in self.coeffs.iter().enumerate() {
                result += *coeff * x.pow(&[i as u64, 0, 0, 0]);
            }

            result
        }
    }

    #[test]
    fn test_extend() {
        let k = 7;

        let coeffs = (0..2usize.pow((k - 1) as u32))
            .map(|i| Fp::from(i as u64))
            .collect::<Vec<Fp>>();

        let poly = UniPoly::new(coeffs);

        let good_curve = GoodCurve::<Fp>::find_k(k);
        let L = prepare_domain(good_curve);

        let s = L[0].iter().step_by(2).map(|x| *x).collect::<Vec<Fp>>();
        let s_prime = L[0]
            .iter()
            .skip(1)
            .step_by(2)
            .map(|x| *x)
            .collect::<Vec<Fp>>();
        let evals = s.iter().map(|s_i| poly.eval(*s_i)).collect::<Vec<Fp>>();

        let (matrices, inverse_matrices) = prepare_matrices(&L);
        let extended_evals = extend(&evals, &L, &matrices, &inverse_matrices, 0);

        // Compute the expected extended evaluations
        let expected_extended_evals = s_prime
            .iter()
            .map(|L_i| poly.eval(*L_i))
            .collect::<Vec<Fp>>();

        assert_eq!(extended_evals, expected_extended_evals);
    }
}
