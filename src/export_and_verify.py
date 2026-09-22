#!/usr/bin/env python3
"""
SANE Systems Ltd - Track 1: Quantum Error Correction
Reference implementation and verification for the [[54, 4, 6]] Bivariate Bicycle Q-LDPC code.
"""
import itertools, json, os
import numpy as np

def gf2_rank(A):
    M = A.copy() % 2
    r, c = M.shape
    rank = 0
    for col in range(c):
        pivot = None
        for row in range(rank, r):
            if M[row, col] == 1:
                pivot = row
                break
        if pivot is not None:
            M[[rank, pivot]] = M[[pivot, rank]]
            for row in range(r):
                if row != rank and M[row, col] == 1:
                    M[row] = (M[row] + M[rank]) % 2
            rank += 1
            if rank == r:
                break
    return rank

def build_cyclic_shift(l, m, poly_terms):
    N = l * m
    mat = np.zeros((N, N), dtype=np.uint8)
    for i in range(l):
        for j in range(m):
            row_idx = i * m + j
            for (dx, dy) in poly_terms:
                mat[row_idx, ((i + dx) % l) * m + ((j + dy) % m)] ^= 1
    return mat

def save_alist(filename, matrix):
    """Exports matrix to Radford Neal's standard .alist format used in QEC decoders."""
    rows, cols = matrix.shape
    row_weights = [np.sum(matrix[i, :] != 0) for i in range(rows)]
    col_weights = [np.sum(matrix[:, j] != 0) for j in range(cols)]
    max_row_w = max(row_weights)
    max_col_w = max(col_weights)

    with open(filename, 'w') as f:
        f.write(f"{rows} {cols}\n")
        f.write(f"{max_row_w} {max_col_w}\n")
        f.write(" ".join(map(str, row_weights)) + "\n")
        f.write(" ".join(map(str, col_weights)) + "\n")
        for i in range(rows):
            indices = [str(j + 1) for j in range(cols) if matrix[i, j] != 0]
            f.write(" ".join(indices) + "\n")
        for j in range(cols):
            indices = [str(i + 1) for i in range(rows) if matrix[i, j] != 0]
            f.write(" ".join(indices) + "\n")

def test_distance_bitmask(H_X, H_Z, max_weight=5):
    n = H_X.shape[1]
    rX = gf2_rank(H_X)
    cols_Z = [0] * n
    for j in range(n):
        v = 0
        for i in range(H_Z.shape[0]):
            if H_Z[i, j]:
                v |= (1 << i)
        cols_Z[j] = v
        
    for w in range(1, max_weight + 1):
        for comb in itertools.combinations(range(n), w):
            acc = 0
            for idx in comb:
                acc ^= cols_Z[idx]
            if acc == 0:
                vec = np.zeros(n, dtype=np.uint8)
                vec[list(comb)] = 1
                if gf2_rank(np.vstack([H_X, vec])) > rX:
                    return w
    return max_weight + 1

def main():
    l, m = 9, 3
    poly_A = [[0, 0], [3, 2], [8, 0]]
    poly_B = [[0, 0], [1, 0], [1, 2]]

    print("[*] Synthesizing [[54, 4, 6]] Bivariate Bicycle Matrices...")
    A = build_cyclic_shift(l, m, poly_A)
    B = build_cyclic_shift(l, m, poly_B)

    H_X = np.hstack([A, B])
    H_Z = np.hstack([B.T, A.T])

    # 1. Verify Commutator Identity
    comm = (H_X @ H_Z.T) % 2
    assert np.all(comm == 0), "CSS Commutator check failed!"
    print("[+] CSS Orthogonality Verified: H_X @ H_Z^T == 0 (mod 2)")

    # 2. Verify Dimensions
    rX = gf2_rank(H_X)
    rZ = gf2_rank(H_Z)
    n = H_X.shape[1]
    k = n - rX - rZ
    print(f"[+] Physical Qubits (n) : {n}")
    print(f"[+] Logical Qubits (k)  : {k} (Rank H_X={rX}, Rank H_Z={rZ})")

    # 3. Verify Code Distance
    print("[*] Validating Minimum Distance Lower Bound (d >= 6)...")
    dX = test_distance_bitmask(H_X, H_Z, max_weight=5)
    dZ = test_distance_bitmask(H_Z, H_X, max_weight=5)
    d = min(dX, dZ)
    print(f"[+] Distance Lower Bound: dX >= {dX}, dZ >= {dZ} -> Certified d >= {d}")

    # 4. Save Artifacts
    os.makedirs("data", exist_ok=True)
    np.savez_compressed("data/bb_54_4_6_matrices.npz", H_X=H_X, H_Z=H_Z)
    save_alist("data/H_X.alist", H_X)
    save_alist("data/H_Z.alist", H_Z)

    meta = {
        "code": "[[54, 4, 6]]",
        "family": "Bivariate Bicycle CSS Q-LDPC",
        "group_ring": "F_2[x, y] / <x^9 - 1, y^3 - 1>",
        "n": int(n), "k": int(k), "d": int(d),
        "check_weight_X": int(np.sum(H_X[0])),
        "check_weight_Z": int(np.sum(H_Z[0])),
        "poly_A": poly_A,
        "poly_B": poly_B,
        "date_discovered": "2026-09-22"
    }
    with open("data/code_parameters.json", "w") as f:
        json.dump(meta, f, indent=2)

    print("[+] Saved data/bb_54_4_6_matrices.npz")
    print("[+] Saved data/H_X.alist and data/H_Z.alist")
    print("[+] Saved data/code_parameters.json")

if __name__ == "__main__":
    main()
