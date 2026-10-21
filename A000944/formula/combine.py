"""Polyhedra U(v,f) = S(v,f)/2 + W_ref(v,f) + W_rr(v,f).

W_rr (rotatory reflections, no fixed cell): if <g> contains no reflection all vertex and face orbits
have even size, so for v or f odd only g = (reflection sigma_h) * (rotation of odd order n >= 3) occur
(symmetry C_nh). With v odd, v = 2 [poles are vertices] + n (a + 2b), a = mirror-vertex orbits,
b = off-mirror orbits; with face poles a >= 1 is odd and a = 1, b = 0 is impossible (all vertices on
the mirror). If v - 2 and v are both prime, only the (v-2)-gonal bipyramid qualifies:
W_rr(v, 2(v-2)) = phi(v-2) / (4(v-2)).  Dually for f odd (prisms).  Otherwise W_rr is taken from
the brute-force statistics (v or f <= 11) or left unknown.

Usage: python combine.py VX FY [nprimes]   (needs sensed_{VX+FY}.json from sensed.py)
"""
import glob
import json
import os
import sys
from fractions import Fraction

import numpy as np

from ms import crt
from refl import wref
from t3 import PRIMES24


def isprime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def phi(n):
    from math import gcd
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


def load_brute_rr():
    stats = {}
    for fn in glob.glob('../brute/stats_*.json'):   # files overlap: merge per cell, do not add
        for key, d in json.load(open(fn)).items():
            stats.setdefault(key, {}).update(d)
    rr = {}
    vmax = 0
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        vmax = max(vmax, v)
        rr[(v, f)] = sum((Fraction(w) for t, w in d.items() if json.loads(t)[0] == 'rr'), Fraction(0))
    return rr, vmax


def wrr(v, f, brute, bvmax):
    """Exact W_rr(v,f) or None if unknown."""
    if v <= bvmax:
        return brute.get((v, f), Fraction(0))
    if f <= bvmax:
        return brute.get((f, v), Fraction(0))
    for a, b in ((v, f), (f, v)):
        if a % 2 == 1 and isprime(a) and isprime(a - 2):
            n = a - 2
            return Fraction(phi(n), 4 * n) if b == 2 * n else Fraction(0)
    return None


def wref_cached(VX, FY, p):
    fn = f'wref_{VX}_{FY}_{p}.npy'
    if os.path.exists(fn):
        return np.load(fn)
    W = wref(VX, FY, p)
    np.save(fn, W)
    return W


def main():
    VX, FY = int(sys.argv[1]), int(sys.argv[2])
    npr = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    S = {tuple(map(int, k.split(','))): v for k, v in json.load(open(f'sensed_{VX + FY}.json')).items()}
    primes = PRIMES24[:npr]
    Ws = [wref_cached(VX, FY, p) for p in primes]
    brute, bvmax = load_brute_rr()
    U = {}
    for v in range(4, VX + 1):
        for f in range(4, FY + 1):
            r = wrr(v, f, brute, bvmax)
            if r is None:
                continue
            res = []
            for W, p in zip(Ws, primes):
                s = S.get((v, f), 0) % p * pow(2, p - 2, p) % p
                rr = r.numerator % p * pow(r.denominator % p, p - 2, p) % p
                res.append((s + int(W[v, f]) + rr) % p)
            u = crt(res, primes)
            assert u == crt(res[:-1], primes[:-1]), ('CRT unstable', v, f)
            assert u >= 0, (v, f, u)
            if u:
                U[(v, f)] = u
    json.dump({f'{v},{f}': u for (v, f), u in sorted(U.items())}, open(f'U_{VX}_{FY}.json', 'w'), indent=0)
    # duality U(v,f) = U(f,v): the half-map grammar is not symmetric in v and f, so this is a real test
    dual_ok = dual_bad = 0
    for (v, f), u in U.items():
        if v < f and f <= VX and v <= FY:
            if wrr(f, v, brute, bvmax) is None:
                continue
            if U.get((f, v), 0) == u:
                dual_ok += 1
            else:
                dual_bad += 1
                print('DUALITY mismatch', v, f, u, U.get((f, v), 0))
    print('duality pairs ok', dual_ok, 'bad', dual_bad)
    vals = [int(l.split()[1]) for l in open('../data/b212438.txt') if l.strip() and not l.startswith('#')]
    it = iter(vals)
    ok = bad = 0
    for n in range(4, 18):
        for k in range(4, 2 * n - 3):
            a = next(it)
            if k <= VX and n <= FY and wrr(k, n, brute, bvmax) is not None:
                got = U.get((k, n), 0)
                if got == a:
                    ok += 1
                else:
                    bad += 1
                    print('A212438 mismatch faces', n, 'vertices', k, got, a)
    print('A212438 cells checked', ok + bad, 'ok', ok, 'bad', bad)
    known = {int(l.split()[0]): int(l.split()[1]) for l in open('../data/b000944.txt') if l.strip() and not l.startswith('#')}
    for n in range(4, VX + 1):
        fs = range(n // 2 + 2, 2 * n - 3)
        if 2 * n - 4 > FY or any(wrr(n, f, brute, bvmax) is None for f in fs):
            continue
        a = sum(U.get((n, f), 0) for f in fs)
        tag = ('ok' if known[n] == a else f'MISMATCH (OEIS {known[n]})') if n in known else 'NEW'
        print(f'A000944({n}) = {a}  {tag}')


if __name__ == '__main__':
    main()
