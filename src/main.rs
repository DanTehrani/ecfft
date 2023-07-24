#![allow(non_snake_case)]
mod curve;
// mod extend;
mod get_2sylow_subgroup;
mod preprocess;
mod utils;

use ff::PrimeField;

// Secp256k1's base field
#[derive(PrimeField)]
#[PrimeFieldModulus = "115792089237316195423570985008687907853269984665640564039457584007908834671663"]
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 5]);

fn main() {
    todo!()
}
