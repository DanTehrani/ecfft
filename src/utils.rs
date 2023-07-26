use halo2curves::ff::PrimeField;

pub fn is_quad_residue<F: PrimeField>(a: F) -> bool {
    let sqrt = a.sqrt();
    sqrt.is_some().into()
}
