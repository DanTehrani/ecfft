use ff::PrimeField;

pub fn is_quad_residue<F: PrimeField>(a: F) -> bool {
    let sqrt = a.sqrt();
    sqrt.is_some().into()
}

#[cfg(test)]
pub(crate) mod tests {
    use ff::PrimeField;

    // Secp256k1's base field
    #[derive(PrimeField)]
    #[PrimeFieldModulus = "115792089237316195423570985008687907853269984665640564039457584007908834671663"]
    #[PrimeFieldGenerator = "3"]
    #[PrimeFieldReprEndianness = "little"]
    pub struct Fp([u64; 5]);
}
