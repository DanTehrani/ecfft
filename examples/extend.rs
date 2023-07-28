use ecfft::{extend, prepare_domain, prepare_matrices, GoodCurve};

type F = halo2curves::bn256::Fq;

fn main() {
    let k = 6;
    let curve = GoodCurve::find_k(k);
    let L = prepare_domain(curve);
    let (matrices, inverse_matrices) = prepare_matrices(&L);

    let evals = (0..(2usize.pow((k - 1) as u32)))
        .map(|i| F::from(i as u64))
        .collect::<Vec<F>>();

    let _extended_evals = extend(&evals, &L, &matrices, &inverse_matrices, 0);
}
