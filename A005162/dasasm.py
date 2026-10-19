"""Count n x n DASASMs (A005162) by a column-sweep transfer matrix on the fundamental triangle.

A DASASM is determined by its entries in T = {(i, j) : i <= j <= n+1-i}. Row i of the full matrix reads
column i of T downward (rows < i, by A = A^T), then row i of T, then column n+1-i of T upward (by antidiagonal
symmetry). So the ASM conditions become a six-vertex model on T:
  * interior cell (i < min(j, n+1-j)): horizontal running sum h and vertical running sum v, both in {0,1};
    entry x with h+x, v+x in {0,1}. (0,0) -> (0,0) or (1,1); (1,1) -> (1,1) or (0,0); (0,1), (1,0) fixed.
  * diagonal cell (j, j), j < n+1-j: the column-j flow v turns into row j; the new row sum r = v + x is free.
  * antidiagonal cell (n+1-j, j), j > n+1-j: row n+1-j ends; x = 1 - h - v is forced and always valid.
  * centre cell (n odd): x = 1 - 2v, always valid.
Sweeping columns left to right, the open rows form a stack (row j opens at column j, closes at column n+1-j),
so the frontier has at most about n/2 bits.

Usage: python dasasm.py NMAX          exact counts for n = 1..NMAX (Python ints; fine to n ~ 40)
"""
import sys


def count(n):
    # f maps a tuple of open-row bits (row 1 first) to a count.
    f = {(): 1}
    for c in range(1, n + 1):
        h = min(c, n + 1 - c)
        # attach the vertical bit v = 0 at the top of column c
        g = {k + (0,): w for k, w in f.items()}
        for r in range(h - 1):  # interior rows 1..h-1 (0-based index r)
            new = {}
            for k, w in g.items():
                hb, v = k[r], k[-1]
                if hb == v:
                    for b in (0, 1):
                        k2 = k[:r] + (b,) + k[r + 1:-1] + (b,)
                        new[k2] = new.get(k2, 0) + w
                else:
                    new[k] = new.get(k, 0) + w
            g = new
        f = {}
        if c < n + 1 - c:      # diagonal: drop v, open row c with a free bit
            for k, w in g.items():
                for b in (0, 1):
                    k2 = k[:-1] + (b,)
                    f[k2] = f.get(k2, 0) + w
        elif c > n + 1 - c:    # antidiagonal: close the top open row and drop v
            for k, w in g.items():
                k2 = k[:-2]
                f[k2] = f.get(k2, 0) + w
        else:                  # centre
            for k, w in g.items():
                k2 = k[:-1]
                f[k2] = f.get(k2, 0) + w
    return f[()]


def odd_formula(m):
    """|DASASM(2m+1)| = prod_{i=0}^m (3i)!/(m+i)!  (Behrend-Fischer-Konvalinka 2017)."""
    from math import factorial
    num = den = 1
    for i in range(m + 1):
        num *= factorial(3 * i)
        den *= factorial(m + i)
    assert num % den == 0
    return num // den


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    for n in range(1, nmax + 1):
        a = count(n)
        chk = ""
        if n % 2:
            chk = "formula ok" if a == odd_formula(n // 2) else "FORMULA MISMATCH"
        print(n, a, chk, flush=True)
