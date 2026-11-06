"""Planar-graph tower mod p (fast version of tower.py; same equations, see its docstring).

Series in the vertex power sums are dense int64 vectors indexed by all partitions of weight <= N.
The 3-connected input is a list of rows (n, vertex_partition, edge_orbits, weight mod p).  Edge
substitutions are expanded in binomial moments: an orbit kind (k, flipped) occurring c times gives
(1 + d_kind)^c = sum_j C(c, j) d_kind^j with d_kind = D^(id/tau)[k] - 1 of weight >= k, so only
j <= (N - n)/k matter.  Rows are compressed once into moment tables; every fixed-point pass then
needs only a few hundred products.

Usage (validation against exact brute-force data):
    python tower_p.py N zk_file.json [...]          -> prints the three sequences (CRT over primes)
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction
from math import comb, gcd

import numpy as np

PRIMES24 = [16777213, 16777199, 16777183, 16777153, 16777141, 16777139, 16777127, 16777121, 16777099,
            16777049, 16777027, 16776989]


def partitions_upto(N):
    out = [()]

    def rec(rem, mx, cur):
        for k in range(min(rem, mx), 0, -1):
            nc = cur + (k,)
            out.append(nc)
            rec(rem - k, k, nc)
    rec(N, N, ())
    return out


class PS:
    """Truncated cycle-index series ring mod p."""

    _cache = {}

    def __init__(self, N, p):
        self.N, self.p = N, p
        key = N
        if key not in PS._cache:
            PS._cache[key] = self._build(N)
        (self.parts, self.idx, self.wt, self.I, self.J, self.K, self.Kstart,
         self.pl, self.d1src, self.d1dst, self.d1mul) = PS._cache[key]
        self.P = len(self.parts)
        # weight caps: outputs of weight <= w are a prefix in K order
        nP = np.searchsorted(self.wt, np.arange(N + 1), side='right')
        self.cap_pairs = np.searchsorted(self.K, nP)          # pairs with K < nP[w]
        self.cap_starts = np.searchsorted(self.Kstart, self.cap_pairs)
        self.diag = self.I == self.J
        self.set_cap(N)

    def set_cap(self, w):
        """Compute products only up to total weight w (entries above are left zero)."""
        self.w = w
        n = int(self.cap_pairs[w])
        ns = int(self.cap_starts[w])
        self._I, self._J, self._diag = self.I[:n], self.J[:n], self.diag[:n]
        self._Ks = self.Kstart[:ns]
        self._Kout = self.K[self._Ks]

    @staticmethod
    def _build(N):
        parts = partitions_upto(N)
        parts.sort(key=lambda l: (sum(l), l))
        idx = {l: i for i, l in enumerate(parts)}
        wt = np.array([sum(l) for l in parts])
        I, J, K = [], [], []
        byw = defaultdict(list)
        for i, l in enumerate(parts):
            byw[sum(l)].append(i)
        for i, a in enumerate(parts):
            wa = sum(a)
            for w in range(0, N - wa + 1):
                for j in byw[w]:
                    if j < i:
                        continue
                    k = idx[tuple(sorted(a + parts[j], reverse=True))]
                    I.append(i)
                    J.append(j)
                    K.append(k)
        I, J, K = map(np.array, (I, J, K))
        order = np.argsort(K, kind='stable')
        I, J, K = I[order], J[order], K[order]
        Kstart = np.flatnonzero(np.r_[True, K[1:] != K[:-1]])
        pl = {}
        for k in range(1, N + 1):
            src, dst = [], []
            for i, l in enumerate(parts):
                if sum(l) * k <= N:
                    src.append(i)
                    dst.append(idx[tuple(x * k for x in l)])
            pl[k] = (np.array(src), np.array(dst))
        d1src, d1dst, d1mul = [], [], []
        for i, l in enumerate(parts):
            c = l.count(1)
            if c:
                l2 = list(l)
                l2.remove(1)
                d1src.append(i)
                d1dst.append(idx[tuple(l2)])
                d1mul.append(c)
        return (parts, idx, wt, I, J, K, Kstart, pl, np.array(d1src), np.array(d1dst), np.array(d1mul))

    def zero(self):
        return np.zeros(self.P, dtype=np.int64)

    def one(self):
        a = self.zero()
        a[0] = 1
        return a

    def mono(self, parts):
        a = self.zero()
        key = tuple(sorted(parts, reverse=True))
        if sum(key) <= self.N:
            a[self.idx[key]] = 1
        return a

    def mul(self, a, b):
        p = self.p
        I, J = self._I, self._J
        out = self.zero()
        if not len(I):
            return out
        prod = (a[I] * b[J]) % p
        prod2 = np.where(self._diag, prod, (prod + (a[J] * b[I]) % p) % p)
        sums = np.add.reduceat(prod2, self._Ks)
        out[self._Kout] = sums % p
        return out

    def pleth(self, a, k):
        src, dst = self.pl[k]
        out = self.zero()
        out[dst] = a[src]
        return out

    def dp1(self, a):
        out = self.zero()
        np.add.at(out, self.d1dst, a[self.d1src] * self.d1mul % self.p)
        return out % self.p

    def inv(self, c):
        return pow(int(c) % self.p, self.p - 2, self.p)

    def exp0(self, f):
        assert f[0] % self.p == 0
        r, term = self.one(), self.one()
        for j in range(1, self.N + 1):
            term = self.mul(term, f) * self.inv(j) % self.p
            if not term.any():
                break
            r = (r + term) % self.p
        return r

    def powr(self, a, e):
        r = self.one()
        for _ in range(e):
            r = self.mul(r, a)
        return r


# ---------------------------------------------------------------- 3-connected rows -> moment tables

def moment_table(rows, N, p):
    """rows: (vpart tuple, {kind: count}, weight mod p).  Returns dict
    vpart -> dict (tuple of (kind, j) sorted) -> coefficient mod p."""
    tab = defaultdict(lambda: defaultdict(int))
    for vp, kinds, w in rows:
        budget = N - sum(vp)
        if budget < 0:
            continue
        # enumerate j-vectors with sum_k j*k <= budget
        items = sorted(kinds.items())
        def rec(i, rem, cur, coef):
            if i == len(items):
                key = tuple(cur)
                tab[vp][key] = (tab[vp][key] + coef) % p
                return
            (k, fl), c = items[i]
            for j in range(0, min(c, rem // k) + 1):
                cc = coef * comb(c, j) % p
                rec(i + 1, rem - j * k, cur + ([((k, fl), j)] if j else []), cc)
        rec(0, budget, [], w % p)
    return tab


def rows_from_zk(zk, p):
    """zk rows (n, vp, eo, w) -> three row lists (R, Hid, Htau), weights mod p."""
    R, Hid, Htau = [], [], []
    for n, vp, eo, w in zk:
        kinds = defaultdict(int)
        for e in eo:
            kinds[tuple(e)] += 1
        R.append((tuple(vp), dict(kinds), w))
        c = kinds.get((1, 0), 0)
        if c:
            l = list(vp)
            l.remove(1)
            l.remove(1)
            k2 = dict(kinds)
            k2[(1, 0)] -= 1
            Hid.append((tuple(l), k2, 2 * c * w % p))
        c = kinds.get((1, 1), 0)
        if c:
            l = list(vp)
            l.remove(2)
            k2 = dict(kinds)
            k2[(1, 1)] -= 1
            Htau.append((tuple(l), k2, 2 * c * w % p))
    return R, Hid, Htau


class Core:
    def __init__(self, ps, R, Hid, Htau):
        self.ps = ps
        self.tabs = [moment_table(x, ps.N, ps.p) for x in (R, Hid, Htau)]

    def evaluate(self, Did, Dtau):
        ps, p = self.ps, self.ps.p
        one = ps.one()
        d = {0: (Did - one) % p, 1: (Dtau - one) % p}
        powcache = {}

        def dpow(kind, j):
            key = (kind, j)
            if key not in powcache:
                if j == 0:
                    powcache[key] = one
                else:
                    base = ps.pleth(d[kind[1]], kind[0])
                    powcache[key] = ps.mul(dpow(kind, j - 1), base)
            return powcache[key]

        prodcache = {}

        def prod(key):
            if key not in prodcache:
                if not key:
                    prodcache[key] = one
                else:
                    prodcache[key] = ps.mul(prod(key[:-1]), dpow(*key[-1]))
            return prodcache[key]

        outs = []
        for tab in self.tabs:
            acc = ps.zero()
            for vp, inner in tab.items():
                m = ps.mono(vp)
                s = ps.zero()
                for key, c in inner.items():
                    if c:
                        s = (s + c * prod(key)) % p
                acc = (acc + ps.mul(m, s)) % p
            outs.append(acc)
        return outs  # R, Hid, Htau


# ---------------------------------------------------------------- networks and blocks

def E_id(ps, Mid):
    p = ps.p
    f = ps.zero()
    for k in range(1, ps.N + 1):
        f = (f + ps.pleth(Mid, k) * ps.inv(k)) % p
    return ps.exp0(f)


def E_tau(ps, Mid, Mtau):
    p = ps.p
    f = ps.zero()
    for k in range(1, ps.N + 1):
        f = (f + ps.pleth(Mtau if k % 2 else Mid, k) * ps.inv(k)) % p
    return ps.exp0(f)


def series_nets(ps, Aid, Atau):
    p = ps.p
    p1, p2 = ps.mono([1]), ps.mono([2])
    Sid = ps.zero()
    t = ps.mul(ps.mul(Aid, Aid), p1)
    while t.any():
        Sid = (Sid + t) % p
        t = ps.mul(ps.mul(t, Aid), p1)
    A2 = ps.pleth(Aid, 2)
    p2A2 = ps.mul(p2, A2)
    Stau = ps.zero()
    t = ps.mul(p1, A2)
    while t.any():
        Stau = (Stau + t) % p
        t = ps.mul(t, p2A2)
    t = ps.mul(Atau, p2A2)
    while t.any():
        Stau = (Stau + t) % p
        t = ps.mul(t, p2A2)
    return Sid, Stau


def solve_networks(ps, core, verbose=False):
    p = ps.p
    one = ps.one()
    Z = ps.zero()
    Sid = Stau = Pid = Ptau = Hid = Htau = Z
    # weight by weight: values below weight w are final once the system is stable at cap w - 1
    for w in range(1, ps.N + 1):
        ps.set_cap(w)
        Sid, Stau, Pid, Ptau, Hid, Htau = _passes(ps, core, (Sid, Stau, Pid, Ptau, Hid, Htau), verbose, w)
    return Sid, Stau, Pid, Ptau, Hid, Htau


def _passes(ps, core, state, verbose, w):
    p = ps.p
    one = ps.one()
    Sid, Stau, Pid, Ptau, Hid, Htau = state
    for it in range(3 * ps.N + 6):
        Did = (one + Sid + Pid + Hid) % p
        Dtau = (one + Stau + Ptau + Htau) % p
        Aid, Atau = (one + Pid + Hid) % p, (one + Ptau + Htau) % p
        Mid, Mtau = (Sid + Hid) % p, (Stau + Htau) % p
        nSid, nStau = series_nets(ps, Aid, Atau)
        nPid = (2 * E_id(ps, Mid) - 2 * one - Mid) % p
        nPtau = (2 * E_tau(ps, Mid, Mtau) - 2 * one - Mtau) % p
        _, nHid, nHtau = core.evaluate(Did, Dtau)
        new = (nSid, nStau, nPid, nPtau, nHid, nHtau)
        same = all((a == b).all() for a, b in zip(new, (Sid, Stau, Pid, Ptau, Hid, Htau)))
        Sid, Stau, Pid, Ptau, Hid, Htau = new
        if same:
            break
    else:
        raise RuntimeError('network fixed point did not converge')
    if verbose:
        print(f'  weight {w}: {it + 1} passes', flush=True)
    return Sid, Stau, Pid, Ptau, Hid, Htau


def blocks(ps, core, verbose=False):
    p = ps.p
    one = ps.one()
    Sid, Stau, Pid, Ptau, Hid, Htau = solve_networks(ps, core, verbose)
    Did = (one + Sid + Pid + Hid) % p
    Dtau = (one + Stau + Ptau + Htau) % p
    Aid, Atau = (one + Pid + Hid) % p, (one + Ptau + Htau) % p
    Mid, Mtau = (Sid + Hid) % p, (Stau + Htau) % p
    p1, p2, p11 = ps.mono([1]), ps.mono([2]), ps.mono([1, 1])
    h = ps.inv(2)
    K2 = (p11 + p2) * h % p
    R, _, _ = core.evaluate(Did, Dtau)
    N = ps.N
    Snode = ps.zero()
    A2 = ps.pleth(Aid, 2)
    for k in range(3, N + 1):
        acc = ps.zero()
        for m in range(1, k + 1):
            if k % m == 0:
                ph = sum(1 for j in range(1, m + 1) if gcd(j, m) == 1)
                acc = (acc + ph * ps.mul(ps.powr(ps.mono([m]), k // m), ps.powr(ps.pleth(Aid, m), k // m))) % p
        if k % 2:
            t = ps.mul(ps.mul(p1, ps.powr(p2, (k - 1) // 2)), ps.mul(Atau, ps.powr(A2, (k - 1) // 2)))
            acc = (acc + k * t) % p
        else:
            t = ps.mul(ps.mul(p11, ps.powr(p2, (k - 2) // 2)), ps.powr(A2, k // 2))
            acc = (acc + (k // 2) * t) % p
            t = ps.mul(ps.powr(p2, k // 2), ps.mul(ps.mul(Atau, Atau), ps.powr(A2, (k - 2) // 2)))
            acc = (acc + (k // 2) * t) % p
        Snode = (Snode + acc * ps.inv(2 * k)) % p
    E2id = (ps.mul(Mid, Mid) + ps.pleth(Mid, 2)) * h % p
    E2tau = (ps.mul(Mtau, Mtau) + ps.pleth(Mid, 2)) * h % p
    Qid = (2 * E_id(ps, Mid) - 2 * one - 2 * Mid - E2id) % p
    Qtau = (2 * E_tau(ps, Mid, Mtau) - 2 * one - 2 * Mtau - E2tau) % p
    Pnode = (ps.mul(p11, Qid) + ps.mul(p2, Qtau)) * h % p

    def pair(Xid, Xtau, Yid, Ytau):
        return (ps.mul(p11, ps.mul(Xid, Yid)) + ps.mul(p2, ps.mul(Xtau, Ytau))) * h % p

    SP = pair(Sid, Stau, Pid, Ptau)
    SR = pair(Sid, Stau, Hid, Htau)
    PR = pair(Pid, Ptau, Hid, Htau)
    RR = (ps.mul((p11 + p2) % p, ps.pleth(Hid, 2)) - ps.mul(p11, ps.mul(Hid, Hid))
          - ps.mul(p2, ps.mul(Htau, Htau))) * ps.inv(4) % p
    B = (K2 + R + Snode + Pnode - SP - SR - PR + RR) % p
    return B, K2


# ---------------------------------------------------------------- univariate stage

def umul(a, b, N, p):
    r = np.zeros(N + 1, dtype=object)
    for i in range(N + 1):
        if a[i]:
            r[i:] += a[i] * b[:N + 1 - i]
    return np.array([int(x) % p for x in r], dtype=np.int64)


def upleth(f, k, N):
    g = np.zeros(N + 1, dtype=np.int64)
    g[0:N + 1:k] = f[:N // k + 1]
    return g


def evalp(ps, Z, f):
    N, p = ps.N, ps.p
    g = {j: upleth(f, j, N) for j in range(1, N + 1)}
    vals = [None] * ps.P
    out = np.zeros(N + 1, dtype=np.int64)
    for i, lam in enumerate(ps.parts):
        if not lam:
            v = np.zeros(N + 1, dtype=np.int64)
            v[0] = 1
        else:
            v = umul(vals[ps.idx[lam[:-1]]], g[lam[-1]], N, p)
        vals[i] = v
        if Z[i]:
            out = (out + Z[i] * v) % p
    return out


def uexp(f, N, p):
    r = [0] * (N + 1)
    r[0] = 1
    df = [i * int(f[i]) % p for i in range(N + 1)]
    for n in range(1, N + 1):
        r[n] = sum(df[k] * r[n - k] for k in range(1, n + 1)) % p * pow(n, p - 2, p) % p
    return np.array(r, dtype=np.int64)


def connected(ps, B):
    N, p = ps.N, ps.p
    Bp = ps.dp1(B)
    C = np.zeros(N + 1, dtype=np.int64)
    C[1] = 1
    for it in range(N + 2):
        bp = evalp(ps, Bp, C)
        s = np.zeros(N + 1, dtype=np.int64)
        for k in range(1, N + 1):
            s = (s + upleth(bp, k, N) * pow(k, p - 2, p)) % p
        e = uexp(s, N, p)
        newC = np.r_[0, e[:N]]
        if (newC == C).all():
            break
        C = newC
    BC = evalp(ps, B, C)
    BpC = evalp(ps, Bp, C)
    Cu = (C + BC - umul(C, BpC, N, p)) % p
    return C, Cu


def euler(Cu, N, p):
    s = np.zeros(N + 1, dtype=np.int64)
    for k in range(1, N + 1):
        s = (s + upleth(Cu, k, N) * pow(k, p - 2, p)) % p
    return uexp(s, N, p)


def run_mod(N, p, rows_fn, verbose=False):
    ps = PS(N, p)
    R, Hid, Htau = rows_fn(p)
    core = Core(ps, R, Hid, Htau)
    B, K2 = blocks(ps, core, verbose)
    one_x = np.zeros(N + 1, dtype=np.int64)
    one_x[1] = 1
    bc = evalp(ps, (B - K2) % p, one_x)
    C, Cu = connected(ps, B)
    G = euler(Cu, N, p)
    return bc, Cu, G


def crt(residues, primes):
    x, m = 0, 1
    for r, q in zip(residues, primes):
        t = (r - x) * pow(m, -1, q) % q
        x += m * t
        m *= q
    return x - m if x > m // 2 else x


def frac_mod(fr, p):
    fr = Fraction(fr)
    return fr.numerator % p * pow(fr.denominator % p, p - 2, p) % p


def main():
    N = int(sys.argv[1])
    files = sys.argv[2:]
    rows = []
    for fn in files:
        rows += [r for r in json.load(open(fn)) if r[0] <= N]

    def rows_fn(p):
        zk = [(n, vp, eo, frac_mod(w, p)) for n, vp, eo, w in rows]
        return rows_from_zk(zk, p)

    nprimes = 3
    res = [run_mod(N, q, rows_fn) for q in PRIMES24[:nprimes]]
    for name, k in (('A021103', 0), ('A003094', 1), ('A005470', 2)):
        vals = [crt([int(r[k][n]) for r in res], PRIMES24[:nprimes]) for n in range(N + 1)]
        print(name, vals)


if __name__ == '__main__':
    main()
