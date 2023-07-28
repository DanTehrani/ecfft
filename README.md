# Elliptic Curve Fast Fourier Transform 

This library implements the EXTEND operation from [ECFFT Part I](https://arxiv.org/pdf/2107.08473.pdf). It includes the following preprocessing algorithms. 

- Finding an elliptic curve over a given odd-sized finite field with a subgroup of order $2^k$.
- Constructing the FFTree based on the 2-isogeny chain of elliptic curves.

## Run tests
```bash
cargo test
```