"""Prototype: count minimal triangle graphs (A000080) by species theory.

Incidence graph I (points X, triples Y): the triple system is minimal iff every Y-node is a cut vertex of I.
Blocks of I are bridges or subdivided 2-connected loopless multigraphs (edges = triples).  So the system is a
point/unit tree, unit = lone triple E3(X) or a 2-connected multigraph block B(X,Y) with Y -> arm point.
B comes from the connected edge-labelled multigraph species C via the block decomposition, evaluated through
plethystic substitution into one-variable series.
"""
import sys
from fractions import Fraction as F
from math import gcd, factorial
from functools import lru_cache

N = int(sys.argv[1]) if len(sys.argv) > 1 else 14


def mul(a, b):
    c = [F(0)] * (N + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(N + 1 - i):
                if b[j]:
                    c[i + j] += ai * b[j]
    return c


def add(*s):
    return [sum(t) for t in zip(*s)]


def scale(a, k):
    return [k * x for x in a]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def at(a, k):  # a(x^k)
    c = [F(0)] * (N + 1)
    for i in range(0, N // k + 1):
        c[i * k] = a[i]
    return c


def exp_(a):  # a[0] == 0
    e = [F(0)] * (N + 1)
    e[0] = F(1)
    for n in range(1, N + 1):
        e[n] = sum(k * a[k] * e[n - k] for k in range(1, n + 1)) / n
    return e


def log_(a):  # a[0] == 1
    l = [F(0)] * (N + 1)
    for n in range(1, N + 1):
        l[n] = (n * a[n] - sum(k * l[k] * a[n - k] for k in range(1, n))) / n
    return l


def inv(a):
    b = [F(0)] * (N + 1)
    b[0] = 1 / a[0]
    for n in range(1, N + 1):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b


def mobius(n):
    r, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            r = -r
        p += 1
    return -r if n > 1 else r


def pexp(a):
    return exp_(add(*[scale(at(a, k), F(1, k)) for k in range(1, N + 1)]))


def plog(a):  # a[0] == 1
    l = log_(a)
    return add(*[scale(at(l, k), F(mobius(k), k)) for k in range(1, N + 1)])


def partitions(n, m=None):
    if m is None:
        m = n
    if n == 0:
        yield []
        return
    for k in range(min(n, m), 0, -1):
        for p in partitions(n - k, k):
            yield [k] + p


def zee(lam):
    z, cnt = 1, {}
    for a in lam:
        cnt[a] = cnt.get(a, 0) + 1
    for a, m in cnt.items():
        z *= a ** m * factorial(m)
    return z, cnt


def pair_orbits(lam):
    o = {}
    for i, a in enumerate(lam):
        if a % 2:
            o[a] = o.get(a, 0) + (a - 1) // 2
        else:
            o[a] = o.get(a, 0) + a // 2 - 1
            o[a // 2] = o.get(a // 2, 0) + 1
        for b in lam[i + 1:]:
            g = gcd(a, b)
            o[a * b // g] = o.get(a * b // g, 0) + g
    return o


PARTS = [(lam,) + zee(lam) + (pair_orbits(lam),) for n in range(N + 1) for lam in partitions(n)]


def graph_series(S, W):
    """Z_G, p1 d/dp1 Z_G, d/dq1 Z_G evaluated at p_k = S(x^k), q_k = W(x^k)."""
    L = add(*[scale(at(W, j), F(1, j)) for j in range(1, N + 1)])  # log of E(W)
    Lat = {d: at(L, d) for d in range(1, N + 1)}
    Sat = {k: at(S, k) for k in range(1, N + 1)}
    G = [F(0)] * (N + 1); Gp = G[:]; Gq = G[:]
    for lam, z, cnt, orb in PARTS:
        t = [F(0)] * (N + 1); t[0] = F(1, z)
        for a in lam:
            t = mul(t, Sat[a])
        if not any(t):
            continue
        ex = [F(0)] * (N + 1)
        for d, m in orb.items():
            if d <= N and m:
                ex = add(ex, scale(Lat[d], m))
        t = mul(t, exp_(ex))
        m1, m2 = cnt.get(1, 0), cnt.get(2, 0)
        G = add(G, t)
        Gp = add(Gp, scale(t, m1))
        Gq = add(Gq, scale(t, m1 * (m1 - 1) // 2 + m2))
    return G, Gp, Gq


def E2(R):
    return scale(add(mul(R, R), at(R, 2)), F(1, 2))


def E3(R):
    return scale(add(mul(mul(R, R), R), scale(mul(R, at(R, 2)), 3), scale(at(R, 3), 2)), F(1, 6))


x = [F(0)] * (N + 1); x[1] = F(1)
R = x[:]; S = x[:]
for it in range(N + 2):
    G, Gp, Gq = graph_series(S, R)
    S = add(S, sub(R, mul(Gp, inv(G))))
    G, Gp, Gq = graph_series(S, R)
    Gi = inv(G)
    Cdot = mul(Gp, Gi)
    CY = mul(Gq, Gi)
    lg = log_(G)
    C = add(*[scale(at(lg, k), F(mobius(k), k)) for k in range(1, N + 1)])
    # B*_X(R,R) = PLog(R/S)
    q = [F(0)] * (N + 1)
    for i in range(N):
        q[i] = R[i + 1]
    s = [F(0)] * (N + 1)
    for i in range(N):
        s[i] = S[i + 1]
    BXs = plog(mul(q, inv(s)))
    BYs = CY
    Bs = add(sub(C, R), mul(R, BXs))
    BX = sub(BXs, mul(R, R))
    BY = sub(BYs, E2(R))
    B = sub(Bs, mul(E2(R), R))
    Ud = add(E2(R), BX, BY)
    U = add(E3(R), B)
    Rn = mul(x, pexp(Ud))
    if Rn == R:
        break
    R = Rn
H = sub(add(R, U), mul(R, Ud))
print("iterations", it)
print("H", [str(c) for c in H])
print("B units", [str(c) for c in B])
