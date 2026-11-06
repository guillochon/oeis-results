"""Planar-graph tower: 3-connected cycle index -> networks -> 2-connected -> connected -> all.

Input: the two-sort cycle index of 3-connected planar graphs, as rows
    [n, vertex_partition, edge_orbits [(k, flipped)], weight]
(brute/zk.py, later the formula version).  Output: A021103 (2-connected), A003094 (connected),
A005470 (all) for n <= N.  Pilot implementation with exact Fractions; series in the vertex
power sums p_1, p_2, ... truncated at total weight N (dict partition -> coefficient).

Networks are Z2-species (pole swap): every network class F carries F^id (pole-fixing cycle index)
and F^tau (pole-swapping, twisted).  D = e + S + P + H.
  S: series of >= 2 networks from A = D - S;  P: set of >= 2 from M = S + H plus at most one e;
  H: edge-rooted 3-connected core, root edge deleted, other edges -> D.
2-connected (dissymmetry on the SPQR tree):
  Z_B = K2 + R + Snode + Pnode - SP - SR - PR + (RR unordered - RR ordered).
Connected: C* = x exp(sum B'(C*)(x^k)/k), C = C* + B(C*) - C* B'(C*).  All: Euler transform.

Usage: python tower.py N zk_file1.json [zk_file2.json ...]
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache

N = 0


def wt(lam):
    return sum(lam)


def merge(a, b):
    return tuple(sorted(a + b, reverse=True))


class Ser(dict):
    """partition(tuple desc) -> Fraction, truncated at weight N"""

    def __add__(self, o):
        r = Ser(self)
        for k, v in o.items():
            r[k] = r.get(k, 0) + v
        return r.clean()

    def __sub__(self, o):
        r = Ser(self)
        for k, v in o.items():
            r[k] = r.get(k, 0) - v
        return r.clean()

    def scale(self, c):
        return Ser({k: v * c for k, v in self.items()}).clean()

    def clean(self):
        for k in [k for k, v in self.items() if v == 0]:
            del self[k]
        return self

    def __mul__(self, o):
        if not isinstance(o, Ser):
            return self.scale(o)
        r = defaultdict(Fraction)
        bw = defaultdict(list)
        for k, v in o.items():
            bw[wt(k)].append((k, v))
        for ka, va in self.items():
            wa = wt(ka)
            for w in range(0, N - wa + 1):
                for kb, vb in bw.get(w, ()):
                    r[merge(ka, kb)] += va * vb
        return Ser(r).clean()

    def pleth(self, k):
        r = Ser()
        for lam, v in self.items():
            if wt(lam) * k <= N:
                key = tuple(x * k for x in lam)
                r[key] = r.get(key, 0) + v
        return r

    def minwt(self):
        return min((wt(k) for k in self), default=N + 1)

    def dp1(self):
        r = Ser()
        for lam, v in self.items():
            c = lam.count(1)
            if c:
                l2 = list(lam)
                l2.remove(1)
                key = tuple(l2)
                r[key] = r.get(key, 0) + v * c
        return r


ONE = None


def mono(parts):
    return Ser({tuple(sorted(parts, reverse=True)): Fraction(1)}) if wt(parts) <= N else Ser()


def powr(s, e):
    r = Ser(ONE)
    for _ in range(e):
        r = r * s
    return r


def exp0(f):
    """exp(f), f without constant term"""
    assert () not in f
    r = Ser(ONE)
    term = Ser(ONE)
    j = 1
    while True:
        term = (term * f).scale(Fraction(1, j))
        if not term:
            break
        r = r + term
        j += 1
    return r


def load_zk(files):
    rows = []
    for fn in files:
        rows += json.load(open(fn))
    out = []
    for n, vp, eo, w in rows:
        if n <= N:
            out.append((n, tuple(vp), tuple(tuple(e) for e in eo), Fraction(w)))
    return out


def subst_edges(eo, Did, Dtau, skip=None):
    """product over edge orbits of D^id[k] (unflipped) or D^tau[k] (flipped); skip one orbit"""
    cnt = defaultdict(int)
    for e in eo:
        cnt[e] += 1
    if skip is not None:
        cnt[skip] -= 1
    r = Ser(ONE)
    for (k, fl), c in cnt.items():
        base = (Dtau if fl else Did).pleth(k)
        r = r * powr(base, c)
    return r


def remove_parts(vp, parts):
    l = list(vp)
    for p in parts:
        l.remove(p)
    return l


def H_series(zk, Did, Dtau):
    Hid, Htau = Ser(), Ser()
    for n, vp, eo, w in zk:
        nfix = sum(1 for e in eo if e == (1, 0))
        nflip = sum(1 for e in eo if e == (1, 1))
        if nfix:
            m = mono(remove_parts(vp, [1, 1]))
            if m:
                Hid = Hid + (m * subst_edges(eo, Did, Dtau, skip=(1, 0))).scale(2 * nfix * w)
        if nflip:
            m = mono(remove_parts(vp, [2]))
            if m:
                Htau = Htau + (m * subst_edges(eo, Did, Dtau, skip=(1, 1))).scale(2 * nflip * w)
    return Hid, Htau


def R_series(zk, Did, Dtau):
    R = Ser()
    for n, vp, eo, w in zk:
        m = mono(vp)
        if m:
            R = R + (m * subst_edges(eo, Did, Dtau)).scale(w)
    return R


def E_id(Mid):
    f = Ser()
    for k in range(1, N + 1):
        f = f + Mid.pleth(k).scale(Fraction(1, k))
    return exp0(f)


def E_tau(Mid, Mtau):
    f = Ser()
    for k in range(1, N + 1):
        f = f + (Mtau if k % 2 else Mid).pleth(k).scale(Fraction(1, k))
    return exp0(f)


def E2(Mid, Mother):
    return (Mother * Mother + Mid.pleth(2)).scale(Fraction(1, 2))


def series_nets(Aid, Atau):
    p1, p2 = mono([1]), mono([2])
    Sid = Ser()
    t = Aid * Aid * p1
    while t:
        Sid = Sid + t
        t = t * Aid * p1
    A2 = Aid.pleth(2)
    Stau = Ser()
    t = p1 * A2                    # k = 2
    while t:
        Stau = Stau + t
        t = t * p2 * A2
    t = Atau * p2 * A2             # k = 3
    while t:
        Stau = Stau + t
        t = t * p2 * A2
    return Sid, Stau


def solve_networks(zk):
    Sid = Stau = Pid = Ptau = Hid = Htau = Ser()
    for it in range(3 * N + 6):          # S and P alternate: about one weight per two passes
        Did = Ser(ONE) + Sid + Pid + Hid
        Dtau = Ser(ONE) + Stau + Ptau + Htau
        Aid, Atau = Ser(ONE) + Pid + Hid, Ser(ONE) + Ptau + Htau
        Mid, Mtau = Sid + Hid, Stau + Htau
        nSid, nStau = series_nets(Aid, Atau)
        nPid = E_id(Mid).scale(2) - Ser(ONE).scale(2) - Mid
        nPtau = E_tau(Mid, Mtau).scale(2) - Ser(ONE).scale(2) - Mtau
        nHid, nHtau = H_series(zk, Did, Dtau)
        new = (nSid, nStau, nPid, nPtau, nHid, nHtau)
        if new == (Sid, Stau, Pid, Ptau, Hid, Htau):
            break
        Sid, Stau, Pid, Ptau, Hid, Htau = new
    else:
        raise RuntimeError('network fixed point did not converge')
    return Sid, Stau, Pid, Ptau, Hid, Htau


def blocks(zk):
    Sid, Stau, Pid, Ptau, Hid, Htau = solve_networks(zk)
    Did = Ser(ONE) + Sid + Pid + Hid
    Dtau = Ser(ONE) + Stau + Ptau + Htau
    Aid, Atau = Ser(ONE) + Pid + Hid, Ser(ONE) + Ptau + Htau
    Mid, Mtau = Sid + Hid, Stau + Htau
    p1, p2 = mono([1]), mono([2])
    p11 = mono([1, 1])
    half = Fraction(1, 2)
    K2 = (p11 + p2).scale(half)
    R = R_series(zk, Did, Dtau)
    # S-nodes: cycles k >= 3, dihedral group, edges -> A
    Snode = Ser()
    A2 = Aid.pleth(2)
    for k in range(3, N + 1):
        acc = Ser()
        for m in range(1, k + 1):
            if k % m == 0:
                from math import gcd
                phi = sum(1 for j in range(1, m + 1) if gcd(j, m) == 1)
                acc = acc + (powr(mono([m]), k // m) * powr(Aid.pleth(m), k // m)).scale(phi)
        if k % 2:
            acc = acc + (p1 * powr(p2, (k - 1) // 2) * Atau * powr(A2, (k - 1) // 2)).scale(k)
        else:
            acc = acc + (p11 * powr(p2, (k - 2) // 2) * powr(A2, k // 2)).scale(k // 2)
            acc = acc + (powr(p2, k // 2) * Atau * Atau * powr(A2, (k - 2) // 2)).scale(k // 2)
        Snode = Snode + acc.scale(Fraction(1, 2 * k))
    # P-nodes: bond, >= 3 edges from M plus at most one e
    Qid = E_id(Mid).scale(2) - Ser(ONE).scale(2) - Mid.scale(2) - E2(Mid, Mid)
    Qtau = E_tau(Mid, Mtau).scale(2) - Ser(ONE).scale(2) - Mtau.scale(2) - E2(Mid, Mtau)
    Pnode = (p11 * Qid + p2 * Qtau).scale(half)

    def pair(Xid, Xtau, Yid, Ytau):
        return (p11 * Xid * Yid + p2 * Xtau * Ytau).scale(half)

    SP = pair(Sid, Stau, Pid, Ptau)
    SR = pair(Sid, Stau, Hid, Htau)
    PR = pair(Pid, Ptau, Hid, Htau)
    RR = ((p11 + p2) * Hid.pleth(2) - p11 * Hid * Hid - p2 * Htau * Htau).scale(Fraction(1, 4))
    B = K2 + R + Snode + Pnode - SP - SR - PR + RR
    return B, K2


def evalp(Z, f):
    """univariate series: Z with p_j -> f(x^j); f list of coefficients (index = power)"""
    out = [Fraction(0)] * (N + 1)
    cache = {}

    def fj(j):
        if j not in cache:
            g = [Fraction(0)] * (N + 1)
            for i, c in enumerate(f):
                if i * j <= N:
                    g[i * j] += c
            cache[j] = g
        return cache[j]

    for lam, v in Z.items():
        g = [Fraction(0)] * (N + 1)
        g[0] = Fraction(1)
        for j in lam:
            g = umul(g, fj(j))
        for i in range(N + 1):
            out[i] += v * g[i]
    return out


def umul(a, b):
    r = [Fraction(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]:
                    r[i + j] += x * b[j]
    return r


def uexp(f):
    assert f[0] == 0
    r = [Fraction(0)] * (N + 1)
    r[0] = Fraction(1)
    # r' = f' r
    df = [i * f[i] for i in range(N + 1)]
    for n in range(1, N + 1):
        r[n] = sum(df[k] * r[n - k] for k in range(1, n + 1)) / n
    return r


def upleth(f, k):
    g = [Fraction(0)] * (N + 1)
    for i, c in enumerate(f):
        if i * k <= N:
            g[i * k] += c
    return g


def connected(B):
    Bp = B.dp1()
    C = [Fraction(0)] * (N + 1)
    C[1] = Fraction(1)
    for it in range(N + 1):
        bp = None
        s = [Fraction(0)] * (N + 1)
        for k in range(1, N + 1):
            bk = upleth(evalp(Bp, C), k)
            for i in range(N + 1):
                s[i] += bk[i] / k
        e = uexp(s)
        newC = [Fraction(0)] + e[:N]
        if newC == C:
            break
        C = newC
    Cr = C
    BC = evalp(B, Cr)
    BpC = evalp(Bp, Cr)
    CB = umul(Cr, BpC)
    Cu = [Cr[i] + BC[i] - CB[i] for i in range(N + 1)]
    return Cr, Cu


def euler(C):
    s = [Fraction(0)] * (N + 1)
    for k in range(1, N + 1):
        ck = upleth(C, k)
        for i in range(N + 1):
            s[i] += ck[i] / k
    return uexp(s)


def main():
    global N, ONE
    N = int(sys.argv[1])
    ONE = Ser({(): Fraction(1)})
    zk = load_zk(sys.argv[2:])
    B, K2 = blocks(zk)
    bc = evalp(B - K2, [Fraction(0), Fraction(1)])
    Cr, Cu = connected(B)
    G = euler(Cu)
    print('A021103', [str(x) for x in bc])
    print('A003094', [str(x) for x in Cu])
    print('A005470', [str(x) for x in G])


if __name__ == '__main__':
    main()
