#![allow(non_snake_case)]
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ecfft::{extend, find_coset_offset, prepare_domain, prepare_matrices, GoodCurve};

type F = ark_secp256k1::Fq;

fn prepare_matrices_bench(c: &mut Criterion) {
    for k in [10] {
        let curve = GoodCurve::<F>::find_k(k);
        let (coset_offset_x, coset_offset_y) =
            find_coset_offset(curve.a, curve.B_sqrt * curve.B_sqrt);
        let L = prepare_domain(
            black_box(curve),
            black_box(coset_offset_x),
            black_box(coset_offset_y),
        );

        let mut group = c.benchmark_group("prepare_matrices");
        group.bench_function(format!("prepare_matrices k = {}", k), |b| {
            b.iter(|| {
                prepare_matrices(black_box(&L));
            })
        });
    }
}

fn prepare_domain_bench(c: &mut Criterion) {
    for k in [10] {
        let curve = GoodCurve::<F>::find_k(k);
        let (coset_offset_x, coset_offset_y) = find_coset_offset(curve.a, curve.B_sqrt);

        let mut group = c.benchmark_group("prepare_domain");
        group.bench_function(format!("prepare_domain k = {}", k), |b| {
            b.iter(|| {
                prepare_domain(
                    black_box(curve),
                    black_box(coset_offset_x),
                    black_box(coset_offset_y),
                );
            })
        });
    }
}

fn extend_bench(c: &mut Criterion) {
    for k in [10] {
        let curve = GoodCurve::<F>::find_k(k);
        let (coset_offset_x, coset_offset_y) = find_coset_offset(curve.a, curve.B_sqrt);
        let L = prepare_domain(curve, coset_offset_x, coset_offset_y);
        let (matrices, inverse_matrices) = prepare_matrices(&L);

        let evals = (0..(2usize.pow((k - 1) as u32)))
            .map(|i| F::from(i as u64))
            .collect::<Vec<F>>();

        let mut group = c.benchmark_group("ecfft");
        group.bench_function(format!("extend k = {}", k), |b| {
            b.iter(|| {
                // Extend the evaluations of a polynomial with degree < 2^(k - 1)
                extend(
                    black_box(&evals),
                    black_box(&matrices),
                    black_box(&inverse_matrices),
                    0,
                );
            })
        });
    }
}

fn set_duration() -> Criterion {
    Criterion::default().sample_size(10)
}

criterion_group! {
    name = benches;
    config = set_duration();
    targets =  extend_bench, prepare_matrices_bench, prepare_domain_bench
}
criterion_main!(benches);
