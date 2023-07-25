use ff::PrimeField;

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

    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
    }

    pub fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            result += *coeff * x.pow(&[i as u64, 0, 0, 0]);
        }

        result
    }
}
