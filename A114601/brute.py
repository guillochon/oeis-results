"""A114601 by brute force: all n x n symmetric positive definite matrices with 2 on the diagonal and
entries in {-1, 0, 1} elsewhere, for small n. Also classifies the *connected* ones (the graph of
nonzero off-diagonal entries is connected) by determinant, which identifies the root lattice they
are a basis of: A_n (det n+1), D_n (det 4), E_6 (det 3), E_7 (det 2), E_8 (det 1).

    python attacks/A114601/brute.py [NMAX]      (default 6; n = 6 takes about a minute)

Method: depth-first over rows. Adding row r to a PD matrix M_k keeps it PD iff the Schur pivot
2 - r^T M_k^{-1} r is positive; the pivot equals det(M_{k+1}) / det(M_k) with both determinants
positive integers, so "pivot >= 1/det(M_k)" is an exact test (up to a safe float margin).
"""
import itertools
import sys
from collections import Counter

import numpy as np

OEIS = [1, 3, 23, 393, 13089, 737595, 58969079]  # a(1..7)


def enumerate_pd(n):
    """Yields (matrix, det) for every PD matrix of size n (as integer numpy arrays)."""
    rows = {k: np.array(list(itertools.product((-1, 0, 1), repeat=k)), dtype=float) for k in range(n)}

    def rec(M, det):
        k = M.shape[0]
        if k == n:
            yield M, det
            return
        inv = np.linalg.inv(M)
        cand = rows[k]
        piv = 2.0 - np.einsum("ij,jk,ik->i", cand, inv, cand)
        for idx in np.nonzero(piv > 0.5 / det)[0]:
            r = cand[idx]
            d = int(round(det * piv[idx]))
            Mn = np.empty((k + 1, k + 1))
            Mn[:k, :k] = M
            Mn[k, :k] = r
            Mn[:k, k] = r
            Mn[k, k] = 2.0
            yield from rec(Mn, d)

    yield from rec(np.array([[2.0]]), 2)


def connected(M):
    n = M.shape[0]
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if j not in seen and M[i, j] != 0:
                seen.add(j)
                stack.append(j)
    return len(seen) == n


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(1, nmax + 1):
        total, conn = 0, Counter()
        for M, det in enumerate_pd(n):
            total += 1
            if connected(M):
                conn[det] += 1
        ok = "OK" if n <= len(OEIS) and total == OEIS[n - 1] else "??"
        print(f"n={n}: total {total} [{ok} vs OEIS]  connected by det: {dict(sorted(conn.items()))}  "
              f"(connected total {sum(conn.values())})")


if __name__ == "__main__":
    main()
