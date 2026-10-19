"""A123234: number of n X n Latin squares up to row and column permutations.

Burnside over S_n x S_n (rows, columns; symbols fixed).

A square fixed by (alpha, beta) has, for each symbol s, its permutation
sigma_s (row -> column of symbol s) satisfying sigma_s alpha = beta sigma_s.
Hence alpha and beta are conjugate, and alpha must be fixed-point-free unless
both are the identity (a fixed row with a moved column would put one symbol
twice in a row).  Then sigma_s ranges over a coset tau*C(alpha) of the
centralizer and the n of them must be pairwise disjoint, so

    Fix(alpha, beta) = n! * D(C(alpha)),

D(G) = number of sharply transitive n-subsets of G.  Moreover, if the cycle
type has parts of different lengths, C(alpha) is a direct product acting on
disjoint supports, and n pairwise disjoint permutations cannot exist on the
smaller support: only types k^m (km = n) contribute.  Therefore

    a(n) = L(n)/(n!)^2 + sum_{k | n, k >= 2} n! * D(C_k wr S_{n/k}) / z^2,
    z = k^m * m!,  L(n) = A002860(n).

D is computed by clique enumeration (D = |G|/n * D1, D1 = sets containing the
identity), except for k = 2 where the fixed squares are 5x10-type arrays and
D(C_2 wr S_m) = 2^(m^2) * M_m / (2m)! with M_m = number of ordered 2m-tuples of
m x m permutation matrices summing to 2J (memoized DP).

Reproduces a(4)..a(8) = 4, 16, 1868, 2420400, 66915816462 (McKay-Meynert-Myrvold).
"""
from math import factorial
from fractions import Fraction
from functools import lru_cache
from itertools import permutations

# A002860: number of Latin squares of order n
L = {
    1: 1, 2: 2, 3: 12, 4: 576, 5: 161280, 6: 812851200, 7: 61479419904000,
    8: 108776032459082956800,
    9: 5524751496156892842531225600,
    10: 9982437658213039871725064756920320000,
    11: 776966836171770144107444346734230682311065600000,
}


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def wreath(k, m):
    """C_k wr S_m acting on km points = centralizer of a permutation of type k^m."""
    n = k * m
    cycles = [list(range(i * k, (i + 1) * k)) for i in range(m)]
    ident = tuple(range(n))
    gens = []
    for c in cycles:
        p = list(ident)
        for j in range(k):
            p[c[j]] = c[(j + 1) % k]
        gens.append(tuple(p))
    for a in range(m):
        for b in range(a + 1, m):
            p = list(ident)
            for x, y in zip(cycles[a], cycles[b]):
                p[x], p[y] = y, x
            gens.append(tuple(p))
    group = {ident}
    frontier = [ident]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                gh = compose(g, h)
                if gh not in group:
                    group.add(gh)
                    new.append(gh)
        frontier = new
    return sorted(group)


def D_clique(G, n):
    """Number of sharply transitive n-subsets of G (pairwise disjoint perms)."""
    ident = tuple(range(n))
    der = [g for g in G if all(x != y for x, y in zip(g, ident))]
    m = len(der)
    adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if all(x != y for x, y in zip(der[i], der[j])):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    bytarget = [0] * n
    for i, g in enumerate(der):
        bytarget[g[0]] |= 1 << i

    def rec(cand, t):  # exactly one element sends point 0 to t, for t = 1..n-1
        if t == n:
            return 1
        c = cand & bytarget[t]
        tot = 0
        while c:
            i = (c & -c).bit_length() - 1
            c &= c - 1
            tot += rec(cand & adj[i], t + 1)
        return tot

    d1 = rec((1 << m) - 1, 1)
    assert (d1 * len(G)) % n == 0
    return d1 * len(G) // n


def M2(m):
    """Ordered 2m-tuples of m x m permutation matrices summing to 2J."""
    perms = list(permutations(range(m)))

    @lru_cache(maxsize=None)
    def f(A, r):
        if r == 0:
            return 1
        tot = 0
        for p in perms:
            B = list(A)
            ok = True
            for i in range(m):
                idx = i * m + p[i]
                if B[idx] == 2:
                    ok = False
                    break
                B[idx] += 1
            if ok:
                tot += f(tuple(B), r - 1)
        return tot

    A = [0] * (m * m)
    for i in range(m):  # first permutation = identity, by column symmetry
        A[i * m + i] = 1
    return factorial(m) * f(tuple(A), 2 * m - 1)


def D_type(k, m):
    n = k * m
    if k == 2:
        num = 2 ** (m * m) * M2(m)
        assert num % factorial(n) == 0
        return num // factorial(n)
    return D_clique(wreath(k, m), n)


def a(n):
    tot = Fraction(L[n], factorial(n) ** 2)
    for k in range(2, n + 1):
        if n % k:
            continue
        m = n // k
        z = k ** m * factorial(m)
        tot += Fraction(factorial(n) * D_type(k, m), z * z)
    assert tot.denominator == 1
    return tot.numerator


if __name__ == "__main__":
    for n in range(4, 12):
        print(n, a(n), flush=True)
