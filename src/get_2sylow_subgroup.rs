use core::panic;

use super::utils::is_quad_residue;
use ff::PrimeField;

fn get_fi_coeffs<F: PrimeField>(xi: F, a: F, b: F) -> ((F, F, F), (F, F, F)) {
    let determinant_sqrt = (xi.square() + a * xi + b).sqrt();
    if determinant_sqrt.is_none().into() {
        panic!("determinant_sqrt {:?} is none", xi);
    } else {
        let f1_deg1_coeff = -F::from(2) * (xi - determinant_sqrt.unwrap());
        let f2_deg1_coeff = -F::from(2) * (xi + determinant_sqrt.unwrap());
        ((F::ONE, f1_deg1_coeff, b), (F::ONE, f2_deg1_coeff, b))
    }
}

fn eval_f1<F: PrimeField>(x: F, a: F, b: F, c: F) -> F {
    a * x.square() + b * x + c
}

// Compute \xi from equation (3) (4) in the paper
fn compute_half_point<F: PrimeField>(x: F, a: F, b: F) -> F {
    // Solve x^2 - 2(\xi - \sqrt{\delta_{\xi}})x + b = 0
    // where \delta_{\xi} = \xi^2 + a\xi + b

    // f_{1, \xi}(x) = x^2 - 2(\xi - \sqrt{\delta_{\xi}})x + b
    let (f1, f2) = get_fi_coeffs(x, a, b);

    // Try solving f_{1, \xi} = 0 for x
    let f1_a = f1.0;
    let f1_b = f1.1;
    let f1_c = f1.2;

    let f2_a = f2.0;
    let f2_b = f2.1;
    let f2_c = f2.2;

    let f1_discriminant = f1_b.square() - F::from(4) * f1_a * f1_c;
    let f2_discriminant = f2_b.square() - F::from(4) * f2_a * f2_c;

    let f1_xi_sqrt_discriminant_sqrt = f1_discriminant.sqrt();
    let f2_xi_sqrt_discriminant_sqrt = f2_discriminant.sqrt();

    if f1_xi_sqrt_discriminant_sqrt.is_none().into()
        && f2_xi_sqrt_discriminant_sqrt.is_none().into()
    {
        panic!("Not roots in f1 nor f2")
    }

    let root = if f1_xi_sqrt_discriminant_sqrt.is_some().into() {
        let root =
            (-f1_b + f1_xi_sqrt_discriminant_sqrt.unwrap()) * (F::from(2) * f1_a).invert().unwrap();
        // Check that the root is correct
        let eval = eval_f1(root, f1_a, f1_b, f1_c);
        debug_assert_eq!(eval, F::ZERO);
        root
    } else {
        let root =
            (-f2_b + f2_xi_sqrt_discriminant_sqrt.unwrap()) * (F::from(2) * f2_a).invert().unwrap();
        // Check that the root is correct
        let eval = eval_f1(root, f2_a, f2_b, f2_c);
        debug_assert_eq!(eval, F::ZERO);
        root
    };

    root
}

// https://www.ams.org/journals/mcom/2005-74-249/S0025-5718-04-01640-0/S0025-5718-04-01640-0.pdf
// Determine the structure of the 2-Sylow cyclic subgroup of
// an elliptic curve
pub fn det_2sylow_cyclic_subgroup<F: PrimeField>(a: F, sqrt_b: F) -> Option<(usize, F)> {
    let b = sqrt_b.square();
    let determinant = a.square() - b.double().double();

    if is_quad_residue(determinant) {
        // If the determinant is a quadratic residue, then the curve
        // has a non-cyclic 2-Sylow subgroup.
        // We don't implement this case.
        None
    } else {
        let mut x = F::ZERO;

        // Check if there is a point that forms a group of order 4
        // We can do this by checking if B is a quadratic residue
        let abscissa_exists = is_quad_residue(b);
        if abscissa_exists {
            // The half point P, 2P = (0, 0)
            x = sqrt_b;
        } else {
            return Some((1, F::ZERO));
        }

        // Loop until x is not a quadratic residue
        let mut k = 2;
        while is_quad_residue(x) {
            let discriminant = (x.square() + a * x + b).sqrt();
            if discriminant.is_none().into() {
                break;
            }

            let x_half = compute_half_point(x, a, b);
            x = x_half;
            k += 1;
        }

        // Return the structure of the 2-Sylow cyclic subgroup
        Some((k, x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fp;
    use ff::Field;

    #[test]
    fn test_det_2sylow_cyclic_subgroup() {
        // Test against the curve with a subgroup of order 2^36 found by https://github.com/andrewmilson/ecfft
        let a = Fp::from_str_vartime(
            "31172306031375832341232376275243462303334845584808513005362718476441963632613",
        )
        .unwrap();
        let sqrt_b = Fp::from_str_vartime(
            "45508371059383884471556188660911097844526467659576498497548207627741160623272",
        )
        .unwrap()
        .sqrt()
        .unwrap();

        let gx = Fp::from_str_vartime(
            "41293412487153066667050767300223451435019201659857889215769525847559135483332",
        )
        .unwrap();

        let (k, x) = det_2sylow_cyclic_subgroup(a, sqrt_b).unwrap();
        assert_eq!(k, 36);
        assert_eq!(x, gx);
    }
}
