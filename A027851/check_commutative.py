"""Independent check of the pi = identity class: its leaves are the commutative semigroups with S^3 != {0}.

Commutative semigroups with S^3 = {0} on n points: generators A (|A| = g), B = S^2 minus {0} (|B| = m),
n = g + m + 1, and the product is any symmetric map A x A -> B u {0} that is onto B.  Up to isomorphism
these are the orbits of S_g x S_m, counted by Burnside:

    N(g, m) = #orbits of all symmetric maps A x A -> B u {0},     onto(g, m) = N(g, m) - N(g, m - 1).

So   #commutative semigroups of order n with S^3 != 0  =  A001426(n) - sum_{g+m+1=n} onto(g, m).

    python check_commutative.py logx_n10.txt [n]
"""
import re, sys
from fractions import Fraction
from math import gcd, factorial

A001426 = [1, 1, 3, 12, 58, 325, 2143, 17291, 221805, 11545843, 3518930337]


def partitions(n, mx=None):
    mx = n if mx is None else mx
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for p in partitions(n - k, k):
            yield (k,) + p


def z(lam):
    r = 1
    for c in set(lam):
        r *= c ** lam.count(c) * factorial(lam.count(c))
    return r


def pair_cycles(lam):
    """cycle lengths of the action on unordered pairs {x, y} (x = y allowed)"""
    out = []
    for i, a in enumerate(lam):
        out += [a] * ((a - 1) // 2)    # pairs inside one cycle, x != y, full orbits
        if a % 2 == 0:
            out.append(a // 2)         # the antipodal pairs
        out.append(a)                  # the pairs {x, x}
        for b in lam[i + 1:]:
            out += [a * b // gcd(a, b)] * gcd(a, b)
    return out


def N(g, m):
    if m < 0:
        return 0
    tot = Fraction(0)
    for al in partitions(g):
        pc = pair_cycles(al)
        for be in partitions(m):
            f = 1
            for L in pc:
                f *= 1 + sum(d for d in be if L % d == 0)
            tot += Fraction(f, z(al) * z(be))
    assert tot.denominator == 1
    return int(tot)


def comm_s3_zero(n):
    return sum(N(g, n - 1 - g) - N(g, n - 2 - g) for g in range(0, n))


if __name__ == "__main__":
    for n in range(1, 11):
        c = comm_s3_zero(n)
        print(n, "commutative with S^3 = 0:", c, " with S^3 != 0:", A001426[n] - c)
    if len(sys.argv) > 1:
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 10
        for line in open(sys.argv[1]):
            m = re.match(r"class \d+ \[([\d,]+)\] .*leaves=(\d+)", line)
            if m and set(m.group(1).split(",")) == {"1"}:
                got = int(m.group(2))
                want = A001426[n] - comm_s3_zero(n)
                print("identity class leaves:", got, "expected:", want, "OK" if got == want else "MISMATCH")
