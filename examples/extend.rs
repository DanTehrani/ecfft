use ecfft::{extend, find_coset_offset, prepare_domain, prepare_matrices, GoodCurve};

type F = halo2curves::bn256::Fq;

fn main() {
    let k = 6;
    let curve = GoodCurve::<F>::find_k(k);
    let (coset_offset_x, coset_offset_y) = find_coset_offset(curve.a, curve.B_sqrt.square());
    let L = prepare_domain(curve, coset_offset_x, coset_offset_y);
    let (matrices, inverse_matrices) = prepare_matrices(&L);

    let evals = (0..(2usize.pow((k - 1) as u32)))
        .map(|i| F::from(i as u64))
        .collect::<Vec<F>>();

    let _extended_evals = extend(&evals, &L, &matrices, &inverse_matrices, 0);
}
