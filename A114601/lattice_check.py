"""Lattice-side check of the A_n and D_n counts behind A114601.

A connected Gram matrix of A114601 is the Gram matrix of a root basis of an irreducible root lattice,
so (connected count for type X) = (number of ordered root bases of X) / |Aut(X)|
                                 = (number of unordered root bases) * n! / |Aut(X)|.
Claimed:
  A_n: root bases <-> oriented spanning trees of K_{n+1}:  #sets = 2^n (n+1)^(n-1),
       |Aut(A_n)| = 2 (n+1)!  (n >= 2)                     ->  c_A(n) = 2^(n-1) (n+1)^(n-2)
  D_n: root bases <-> connected unicyclic signed multigraphs on the n coordinates whose cycle is
       "unbalanced" (the n roots are independent), and every such configuration generates all of
       D_n (index 1):  #sets = U_n 4^n / 2 + (n-1) n^(n-2) 4^(n-1),   U_n = A057500(n),
       |Aut(D_n)| = 2^n n! (n >= 5; 3 times that for n = 4)  ->  c_D(n) = #sets / 2^n (/3 for n = 4)

Here: enumerate all n-subsets of roots directly (A_n: n <= 5, D_n: n <= 5), keep those with Gram
determinant det(X) (a basis of the whole lattice; a sublattice of index k has det(X) k^2), and
compare with the claimed #sets. Also counts the full-rank subsets of index > 1, to confirm that for
D_n they never occur among connected unicyclic configurations... they are reported separately.

    python attacks/A114601/lattice_check.py
"""
import itertools
from math import comb, factorial

import numpy as np


def roots_A(n):
    """Roots e_i - e_j of A_n in R^(n+1)."""
    out = []
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                v = [0] * (n + 1)
                v[i], v[j] = 1, -1
                out.append(v)
    return np.array(out)


def roots_D(n):
    """Roots +-e_i +- e_j of D_n in R^n."""
    out = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * n
                v[i], v[j] = si, sj
                out.append(v)
    return np.array(out)


def unicyclic(n):
    """A057500: labelled connected unicyclic simple graphs on n vertices."""
    return sum(factorial(n) // factorial(n - k) * n ** (n - k) // (2 * n) for k in range(3, n + 1)) if n >= 3 else 0


def count_bases(R, n, det_full):
    """Unordered n-subsets of the root list R (up to sign of each root) that are independent;
    split by Gram determinant."""
    # one representative per +-pair: keep roots whose first nonzero coordinate is positive
    reps = [r for r in R if r[np.nonzero(r)[0][0]] > 0]
    reps = np.array(reps)
    dets = {}
    for S in itertools.combinations(range(len(reps)), n):
        V = reps[list(S)]
        d = int(round(np.linalg.det(V @ V.T)))
        if d:
            dets[d] = dets.get(d, 0) + 1
    # each unordered set of roots (with signs) = subset of representatives times 2^n sign choices
    return {d: c * 2 ** n for d, c in dets.items()}, dets.get(det_full, 0) * 2 ** n


def main():
    print("A_n: bases of the whole lattice (det n+1) vs 2^n (n+1)^(n-1)")
    for n in range(2, 6):
        by_det, full = count_bases(roots_A(n), n, n + 1)
        claim = 2 ** n * (n + 1) ** (n - 1)
        cA = full * factorial(n) // (2 * factorial(n + 1))
        print(f"  n={n}: {full} vs {claim} {'OK' if full == claim else 'MISMATCH'};  c_A = {cA};  "
              f"independent sets by Gram det: {dict(sorted(by_det.items()))}")
    print("D_n: bases of the whole lattice (det 4) vs U_n 4^n/2 + (n-1) n^(n-2) 4^(n-1)")
    for n in range(3, 6):
        by_det, full = count_bases(roots_D(n), n, 4)
        claim = unicyclic(n) * 4 ** n // 2 + (n - 1) * n ** (n - 2) * 4 ** (n - 1)
        aut = 2 ** n * factorial(n) * (3 if n == 4 else 1)
        cD = full * factorial(n) // aut
        print(f"  n={n}: {full} vs {claim} {'OK' if full == claim else 'MISMATCH'};  c_D = {cD};  "
              f"independent sets by Gram det: {dict(sorted(by_det.items()))}")


if __name__ == "__main__":
    main()
