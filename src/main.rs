#![allow(non_snake_case)]
mod get_2sylow_subgroup;
mod utils;

use ff::{Field, PrimeField};
use get_2sylow_subgroup::det_2sylow_cyclic_subgroup;

// Secp256k1's base field
#[derive(PrimeField)]
#[PrimeFieldModulus = "115792089237316195423570985008687907853269984665640564039457584007908834671663"]
#[PrimeFieldGenerator = "3"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 5]);

fn main() {
    type F = Fp;

    let target_k = 10;
    // Find a curve with a 2-Sylow subgroup of order 10
    let mut rng = rand::thread_rng();
    loop {
        let a = F::random(&mut rng);
        let sqrt_b = F::random(&mut rng);
        match det_2sylow_cyclic_subgroup(a, sqrt_b) {
            Some((k, x)) => {
                if k == target_k {
                    println!("a = {:?}", a);
                    println!("sqrt_b = {:?}", sqrt_b);
                    println!("k = {}", k);
                    println!("x = {:?}", x);
                    break;
                }
            }
            None => continue,
        }
    }
}
