"""Brute-force check of a Burnside decomposition of unlabeled posets.

a(n) = (1/n!) sum over cycle types lam of |C_lam| fix(lam).
Claims tested:
  (1) if the cycle lengths of g are pairwise coprime (fixed points count as length 1),
      fix(g) = P(#cycles), the labeled poset count on the quotient;
  (2) for g = two 2-cycles, fix = P(n-2) + 4 K(n-2), K(m) = #labeled posets on [m] with 0 covered by 1.
"""
import sys
from math import factorial, gcd
from itertools import combinations
from collections import Counter

def posets(n):
    """All labeled posets on range(n) as tuples of down-set bitmasks (strict)."""
    res = [()]
    for k in range(n):
        new = []
        for p in res:
            up = [0] * k
            for x in range(k):
                for y in range(k):
                    if p[y] >> x & 1:
                        up[x] |= 1 << y
            for D in range(1 << k):
                # D must be a down-set
                if any(D >> x & 1 and p[x] & ~D for x in range(k)):
                    continue
                for U in range(1 << k):
                    if U & D:
                        continue
                    if any(U >> x & 1 and up[x] & ~U for x in range(k)):
                        continue
                    # every d in D below every u in U, i.e. D subset of down(u) for u in U
                    if any(U >> u & 1 and D & ~p[u] for u in range(k)):
                        continue
                    q = list(p)
                    for u in range(k):
                        if U >> u & 1:
                            q[u] |= 1 << k
                    new.append(tuple(q) + (D,))
        res = new
    return res

def partitions(n, m=None):
    if m is None:
        m = n
    if n == 0:
        yield []
        return
    for k in range(min(n, m), 0, -1):
        for r in partitions(n - k, k):
            yield [k] + r

def perm_of(lam):
    g, s = [], 0
    for c in lam:
        g += [s + (i + 1) % c for i in range(c)]
        s += c
    return g

def class_size(lam, n):
    z = 1
    for c, m in Counter(lam).items():
        z *= c ** m * factorial(m)
    return factorial(n) // z

def invariant(p, g):
    n = len(p)
    for y in range(n):
        for x in range(n):
            if (p[y] >> x & 1) != (p[g[y]] >> g[x] & 1):
                return False
    return True

def covers(p, a, b):
    """a covered by b in p."""
    if not p[b] >> a & 1:
        return False
    return not any(p[b] >> z & 1 and p[z] >> a & 1 for z in range(len(p)))

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
P = {m: posets(m) for m in range(N + 1)}
cnt = {m: len(P[m]) for m in P}
K = {m: sum(covers(p, 0, 1) for p in P[m]) if m >= 2 else 0 for m in P}
print("P:", cnt)
print("K:", K)
for n in range(1, N + 1):
    tot = 0
    for lam in partitions(n):
        g = perm_of(lam)
        fx = sum(invariant(p, g) for p in P[n])
        tot += class_size(lam, n) * fx
        coprime = all(gcd(a, b) == 1 for a, b in combinations(lam, 2))
        note = ""
        if coprime:
            note = "coprime: P(%d)=%d %s" % (len(lam), cnt[len(lam)], "OK" if fx == cnt[len(lam)] else "FAIL")
        elif sorted(lam) == sorted([2, 2] + [1] * (n - 4)):
            pred = cnt[n - 2] + 4 * K[n - 2]
            note = "2+2: P+4K=%d %s" % (pred, "OK" if fx == pred else "FAIL")
        print(n, lam, fx, note)
    print("a(%d) =" % n, tot // factorial(n), tot % factorial(n))
