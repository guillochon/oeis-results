"""A123234 independent check: Burnside's lemma summed over every pair of cycle types.

Reads the output of fixcount.c (lines "lambda;mu;N", N = Fix(alpha, beta)/n! for one pair of those
types) and returns

    a(n) = L(n)/(n!)^2 + n! * sum_{(lambda, mu) != (1^n, 1^n)} N(lambda, mu) / (z_lambda z_mu),

since the class of type lambda has n!/z_lambda elements.  It also reports which pairs fix anything.
A pair reported as "skipped" (only (2^m, 2^m) at n = 10) is filled in from rc_classes.D_type(2, m),
the matrix-tuple formula, which the direct search confirms for m = 2, 3, 4.

    python burnside_check.py logs/fixcount_9.txt logs/fixcount_10_fpf.txt
"""
import sys
from collections import Counter
from fractions import Fraction
from math import factorial

from rc_classes import L, a, D_type


def z(parts):
    out = 1
    for k, m in Counter(parts).items():
        out *= k ** m * factorial(m)
    return out


def check(path):
    lines = [ln.split(";") for ln in open(path).read().split("\n") if ln.strip()]
    first = list(map(int, lines[0][0].split()))
    n = sum(first)
    tot = Fraction(L[n], factorial(n) ** 2)
    nonzero = []
    for la, mu, N in lines:
        la, mu = list(map(int, la.split())), list(map(int, mu.split()))
        if N == "skipped":
            assert set(la) == set(mu) == {2}
            N = D_type(2, n // 2)
            print(f"  ({la} x {mu} not searched; N = {N} from the formula)")
        N = int(N)
        if N:
            nonzero.append((la, mu, N))
            tot += Fraction(factorial(n) * N, z(la) * z(mu))
    assert tot.denominator == 1
    return n, len(lines), nonzero, tot.numerator


if __name__ == "__main__":
    for path in sys.argv[1:]:
        n, pairs, nonzero, val = check(path)
        print(f"n = {n}: {pairs} pairs searched; pairs with fixed squares:")
        for la, mu, N in nonzero:
            print(f"    {la} x {mu}: N = {N}")
        print(f"  Burnside total {val}; rc_classes.a(n) = {a(n)}; {'MATCH' if val == a(n) else 'MISMATCH'}")
