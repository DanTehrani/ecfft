use ecfft::GoodCurve;
use halo2curves::bn256::Fq;

fn main() {
    let k = 6;
    let curve = GoodCurve::<Fq>::find_k(k);

    println!(
        "Found curve: y^2 = x(x^2 + {:?}x + {:?})
with subgroup of order 2^{} with generator:
({:?}, {:?})",
        curve.a,
        curve.B_sqrt.square(),
        curve.k,
        curve.gx,
        curve.gy
    );
}
