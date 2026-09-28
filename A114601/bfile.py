"""Write b114601.txt (n = 1..NMAX) from the proved formula, and check the asymptotic constant.

    python attacks/A114601/bfile.py [NMAX]
"""
import sys
from math import e, exp, gamma, lgamma, log
from pathlib import Path

from formula import C_E, GRAM8, OEIS, a_from_c, connected

HERE = Path(__file__).resolve().parent
GRAM9 = 63545326593  # a(9) by direct enumeration (logs/gram9.log)


def asym_constant():
    """a(n) ~ K n! (2e)^n n^(-7/8). See proofs.md section 5."""
    x = 1 / (2 * e)  # the singularity; T(2x) = 1 there
    reg = (e / 4 - 0.5 + 1 / (4 * e))  # C_A(x0)
    reg += -0.25 - x**2 / 2 - 16 * x**3 / 6 - 208 * x**4 / 24  # C_D(x0) minus its -log(1-t)/4 part
    reg += sum(c * x**k / gamma(k + 1) for k, c in C_E.items())
    return 2 ** (-1 / 8) * exp(reg) / gamma(1 / 8)


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    a = a_from_c([0] + [connected(k) for k in range(1, nmax + 1)], nmax)
    assert a[1:8] == OEIS and a[8] == GRAM8 and a[9] == GRAM9
    (HERE / "b114601.txt").write_text("".join(f"{n} {a[n]}\n" for n in range(1, nmax + 1)))
    K = asym_constant()
    print(f"wrote b114601.txt, n = 1..{nmax}; K = {K:.10f}")
    for n in (10, 25, 50, 100, 200):
        if n <= nmax:
            r = exp(log(a[n]) - lgamma(n + 1) - n * log(2 * e) + 7 / 8 * log(n))
            print(f"  n={n}: a(n) / (n! (2e)^n n^(-7/8)) = {r:.8f}   ratio to K {r / K:.6f}")


if __name__ == "__main__":
    main()
