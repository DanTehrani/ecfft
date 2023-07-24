use crate::curve::{ec_add, ec_mul, AffinePoint, GoodCurve};
use ff::PrimeField;
use num_bigint::BigUint;

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
    let curves = isogeny_chain(good_curve);

    let g = AffinePoint::new(good_curve.gx, good_curve.gy);
    // Generate the group of order 2^k and collect the x-coordinates of the points
    // TODO: Understand why this needs to be k - 1
    let G = (0..2u32.pow((good_curve.k - 1) as u32))
        .map(|i| {
            let s = BigUint::from(i);
            ec_mul::<F, GoodCurve<F>>(&g, &s, &good_curve)
        })
        .collect::<Vec<AffinePoint<F>>>();

    let coset_offset = ec_mul::<F, GoodCurve<F>>(&g, &BigUint::from(100000u32), &good_curve);

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
    use crate::Fp;

    #[test]
    fn test_iso_chain() {
        let good_curve = GoodCurve::<Fp>::find_k(6);
        let curves = isogeny_chain(good_curve);

        // The generator of the group of order 4 should be as expected
        let order_4_gx = curves[curves.len() - 1];
        let B_sqrt = order_4_gx.B_sqrt;
        assert!(order_4_gx.gx == B_sqrt || order_4_gx.gx == -B_sqrt);
    }

    #[test]
    fn test_prepare_domain() {
        let good_curve = GoodCurve::<Fp>::find_k(6);
        let domain = prepare_domain(good_curve);
        for i in 0..domain.len() {
            assert_eq!(domain[i].len(), 2usize.pow((good_curve.k - i - 1) as u32));
        }
    }
}
