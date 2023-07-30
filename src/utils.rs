use halo2curves::ff::PrimeField;

pub fn is_quad_residue<F: PrimeField>(a: F) -> bool {
    let sqrt = a.sqrt();
    sqrt.is_some().into()
}

pub fn find_coset_offset<F: PrimeField>(a: F, b: F) -> (F, F) {
    let mut rng = rand::thread_rng();
    let mut coset_offset_x = F::ZERO;
    let mut coset_offset_y = F::ZERO;

    loop {
        coset_offset_x = F::random(&mut rng);
        let y = (coset_offset_x * (coset_offset_x.square() + a * coset_offset_x + b)).sqrt();

        if y.is_some().into() {
            coset_offset_y = y.unwrap();
            return (coset_offset_x, coset_offset_y);
        }
    }
}
