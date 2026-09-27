"""Turn the output of c/skew.c into A109457, A109458, A109459 and check them.

Usage: python sequences.py skew_output.txt [NMAX]

Structure theorem (see README.md): a Krom function is either empty or is
(k forced variables) + (a partition of the remaining m variables into r
blocks of equivalent literals) + (a skew poset on the r blocks).

  A109457(n) = 1 + sum_k C(n,k) 2^k L(n-k),   L(m) = sum_r S2(m,r) 2^(m-r) T(r)
  A109459(n) = 1 + sum_{m<=n} U(m)
  A109458(n) = 1 + sum_{m<=n} (n-m+1) U'(m)

T(r)  = labelled skew posets on r labelled pairs      (line "S r ... Tall ...")
U(x)  = prod_m (1-x^m)^(-cU(m))   Euler transform of connected counts
U'(x) = prod_m (1-x^m)^(-cUp(m))
"""
import sys
from math import comb, factorial

KNOWN = {
    "A109457": [2, 4, 16, 166, 4170, 224716, 24445368, 5167757614, 2061662323954],
    "A109458": [2, 4, 12, 48, 308, 3028, 49490, 1350894, 62154403],
    "A109459": [2, 3, 6, 14, 45, 196, 1360, 15631],
}


def stirling2(m, r):
    return sum((-1) ** (r - j) * comb(r, j) * j ** m for j in range(r + 1)) // factorial(r)


def euler_transform(c, nmax):
    """U with sum U(m) x^m = prod_m (1-x^m)^(-c[m]), c indexed from 1."""
    u = [0] * (nmax + 1)
    u[0] = 1
    for m in range(1, nmax + 1):
        # multiply by (1-x^m)^(-c[m]) = sum_j C(c+j-1, j) x^(mj)
        cm = c.get(m, 0)
        if cm == 0:
            continue
        new = [0] * (nmax + 1)
        for i in range(nmax + 1):
            if u[i] == 0:
                continue
            j = 0
            while i + m * j <= nmax:
                new[i + m * j] += u[i] * comb(cm + j - 1, j)
                j += 1
        u = new
    return u


def exp_formula(tc, nmax):
    """Labelled: T(r) from connected labelled counts tc[r] (exponential formula)."""
    T = [1] + [0] * nmax
    for n in range(1, nmax + 1):
        # T(n) = sum_{k=1..n} C(n-1,k-1) tc(k) T(n-k)
        T[n] = sum(comb(n - 1, k - 1) * tc.get(k, 0) * T[n - k] for k in range(1, n + 1))
    return T


def main():
    path = sys.argv[1]
    T, Tconn, S, cU, cUp = {0: 1}, {}, {0: 1}, {}, {}
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "S":
            r = int(p[1]); S[r] = int(p[2]); T[r] = int(p[4]); Tconn[r] = int(p[5])
        elif p[0] == "C":
            m = int(p[1]); cU[m] = int(p[2]); cUp[m] = int(p[3])
    rmax = max(T)
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else rmax

    print("skew posets up to iso S(r):", [S[r] for r in sorted(S)])
    print("labelled T(r):            ", [T[r] for r in sorted(T)])
    Texp = exp_formula(Tconn, rmax)
    ok = all(Texp[r] == T[r] for r in range(rmax + 1))
    print("exponential-formula check on T(r) from connected counts:", "OK" if ok else "FAIL " + str(Texp))

    L = {m: sum(stirling2(m, r) * 2 ** (m - r) * T[r] for r in range(0, m + 1)) for m in range(nmax + 1)}
    A109457 = [1 + sum(comb(n, k) * 2 ** k * L[n - k] for k in range(n + 1)) for n in range(nmax + 1)]
    U = euler_transform(cU, nmax)
    Up = euler_transform(cUp, nmax)
    A109459 = [1 + sum(U[m] for m in range(n + 1)) for n in range(nmax + 1)]
    A109458 = [1 + sum((n - m + 1) * Up[m] for m in range(n + 1)) for n in range(nmax + 1)]

    print("U(m) :", U)
    print("U'(m):", Up)
    status = True
    for name, seq in (("A109457", A109457), ("A109458", A109458), ("A109459", A109459)):
        known = KNOWN[name]
        n_ok = min(len(known), len(seq))
        match = seq[:n_ok] == known[:n_ok]
        status &= match
        print(f"{name}: {seq}")
        print(f"   known {known}  ->  {'MATCH' if match else 'MISMATCH'} on n<={n_ok-1}"
              + (f"; NEW: " + ", ".join(f"a({n})={seq[n]}" for n in range(n_ok, len(seq))) if len(seq) > n_ok else ""))
    print("ALL KNOWN TERMS REPRODUCED" if status else "*** SOME KNOWN TERMS DO NOT MATCH ***")
    if "--bfiles" in sys.argv:   # write b-files next to this script (b109457.txt etc.)
        from pathlib import Path
        for name, seq in (("A109457", A109457), ("A109458", A109458), ("A109459", A109459)):
            out = Path(__file__).resolve().parent / f"b{name[1:]}.txt"
            out.write_text("".join(f"{n} {v}\n" for n, v in enumerate(seq)))
            print("wrote", out.name, "n = 0..", len(seq) - 1)


if __name__ == "__main__":
    main()
