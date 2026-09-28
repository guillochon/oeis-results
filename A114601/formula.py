"""A114601 from the root-lattice formula.

a(n) = n! [x^n] exp( sum_{k>=1} c(k) x^k / k! ),   c(k) = c_A(k) + c_D(k) + c_E(k), where
  c_A(1) = 1,  c_A(k) = 2^(k-1) (k+1)^(k-2)                        (k >= 2; A_k)
  c_D(k) = 2^(k-1) U(k) + (k-1) k^(k-2) 2^(k-2)                    (k >= 5; D_k)
  c_D(4) = (2^3 U(4) + 3 * 4^2 * 2^2) / 3 = 104                    (D_4: extra triality symmetry)
  (c_D(3) is not added: D_3 = A_3)
  c_E(6), c_E(7), c_E(8): E_6, E_7, E_8 (constants; c_E(k) = 0 otherwise)
U(k) = A057500(k) = number of labelled connected unicyclic graphs on k vertices.

    python attacks/A114601/formula.py [NMAX]
"""
import sys
from fractions import Fraction
from math import comb, factorial

OEIS = [1, 3, 23, 393, 13089, 737595, 58969079]
GRAM8 = 4925177553  # a(8) by direct enumeration (rust gram 8)

# E-type terms: E_6, E_7 from brute.py (connected, det 3 at n = 6, det 2 at n = 7); all three also
# from the root-basis counts of rust/ (route 2), and E_8 from both routes (logs/roots_e8.log, gram8.log).
C_E = {6: 355104, 7: 43733504, 8: 4053973888}


def unicyclic(k):
    return sum(factorial(k) // factorial(k - j) * k ** (k - j) // (2 * k) for j in range(3, k + 1)) if k >= 3 else 0


def c_A(k):
    return 1 if k == 1 else 2 ** (k - 1) * (k + 1) ** (k - 2)


def c_D(k):
    if k < 4:
        return 0
    v = 2 ** (k - 1) * unicyclic(k) + (k - 1) * k ** (k - 2) * 2 ** (k - 2)
    return v // 3 if k == 4 else v


def connected(k):
    return c_A(k) + c_D(k) + C_E.get(k, 0)


def a_from_c(c, nmax):
    """Labelled exponential formula: a(n) = sum_k C(n-1, k-1) c(k) a(n-k), a(0) = 1."""
    a = [1]
    for n in range(1, nmax + 1):
        a.append(sum(comb(n - 1, k - 1) * c[k] * a[n - k] for k in range(1, n + 1)))
    return a


def c_from_a(a):
    """Inverse: connected counts from the totals."""
    c = [0]
    for n in range(1, len(a)):
        c.append(a[n] - sum(comb(n - 1, k - 1) * c[k] * a[n - k] for k in range(1, n)))
    return c


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    c_data = c_from_a([1] + OEIS)
    print("connected counts from the OEIS data vs formula (A + D + E):")
    for k in range(1, len(OEIS) + 1):
        f = connected(k)
        note = "OK" if f == c_data[k] else f"differs by {c_data[k] - f} (= E-type term not yet included?)"
        print(f"  k={k}: data {c_data[k]}, formula {f}  (A {c_A(k)}, D {c_D(k)}, E {C_E.get(k, 0)})  {note}")
    if 7 not in C_E:
        print(f"  => E_7 term implied by the OEIS a(7): {c_data[7] - c_A(7) - c_D(7)}")
    c = [0] + [connected(k) for k in range(1, nmax + 1)]
    a = a_from_c(c, nmax)
    print("\na(n) from the formula" + ("" if 8 in C_E else " (n >= 8 lacks the E_8 term; n >= 7 lacks E_7 unless filled in)") + ":")
    for n in range(1, nmax + 1):
        mark = ("OK" if a[n] == OEIS[n - 1] else "MISMATCH") if n <= len(OEIS) else "new"
        print(f"  a({n}) = {a[n]}  {mark}")


if __name__ == "__main__":
    main()
