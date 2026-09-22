# [[54, 4, 6]] Sparse Bivariate Bicycle Quantum LDPC Code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22904621.svg)](https://doi.org/10.5281/zenodo.22904621)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Algebraic generators, parity-check matrices, and verification engine for a **[[54, 4, 6]] Bivariate Bicycle (BB) Quantum Low-Density Parity-Check (Q-LDPC) code**.

## Parameters
* **Physical Qubits (n):** 54
* **Logical Qubits (k):** 4
* **Code Distance (d):** 6
* **Logical Rate (k/n):** ~7.41% (~3x higher than a planar surface code at d=6)
* **Check Weight:** Uniform weight-6 checks (w_X = 6, w_Z = 6)
* **Group Ring:** F_2[x, y] / <x^9 - 1, y^3 - 1>

## Polynomial Generators
* A(x, y) = 1 + x^3 y^2 + x^8
* B(x, y) = 1 + x + x y^2

CSS parity check matrices:
H_X = [A | B],  H_Z = [B^T | A^T]
Orthogonality holds identically: H_X @ H_Z^T == 0 (mod 2)

## Verification
Run: `python3 src/export_and_verify.py`
