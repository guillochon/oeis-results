"""Polyhedra U(v,f) = S/2 + W_ref + W_rr on the box v <= VX, f <= FY, every symmetry type included.

Needs sensed_{VX+FY}.json (sensed.py), wref_VX_FY_p.npy (run_wref.py) and wrr_VX_FY_p.npy (run_wrr.py)
for the first nprimes primes of t3.PRIMES24.  Checks: CRT stability, U(v,f) = U(f,v), all of A212438,
A000944(n) for every known n, A002840 (by edges) and A006445 (achiral, by edges) where the box covers
them (using duality to reach cells outside it).  Usage: python combine_full.py VX FY [nprimes]
"""
import json
import sys

import numpy as np

from ms import crt
from t3 import PRIMES24


def bfile(name):
    return {int(l.split()[0]): int(l.split()[1]) for l in open(f'../data/{name}') if l.strip() and not l.startswith('#')}


def main():
    VX, FY = int(sys.argv[1]), int(sys.argv[2])
    npr = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    primes = PRIMES24[:npr]
    S = {tuple(map(int, k.split(','))): v for k, v in json.load(open(f'sensed_{VX + FY}.json')).items()}
    Wr = [np.load(f'wref_{VX}_{FY}_{p}.npy') for p in primes]
    Wq = [np.load(f'wrr_{VX}_{FY}_{p}.npy') for p in primes]
    U = {}
    for v in range(4, VX + 1):
        for f in range(4, FY + 1):
            res = []
            for a, b, p in zip(Wr, Wq, primes):
                s = S.get((v, f), 0) % p * pow(2, p - 2, p) % p
                res.append((s + int(a[v, f]) + int(b[v, f])) % p)
            u = crt(res, primes)
            assert u == crt(res[:-1], primes[:-1]) and u >= 0, ('CRT unstable / negative', v, f, u)
            if u:
                U[(v, f)] = u
    json.dump({f'{v},{f}': u for (v, f), u in sorted(U.items())}, open(f'Ufull_{VX}_{FY}.json', 'w'), indent=0)

    def u(v, f):
        if (v, f) in U or (v <= VX and f <= FY):
            return U.get((v, f), 0)
        if f <= VX and v <= FY:
            return U.get((f, v), 0)
        return None

    dual = sum(1 for (v, f) in U if f <= VX and v <= FY and U.get((f, v), 0) != U[(v, f)])
    npairs = sum(1 for (v, f) in U if v < f and f <= VX and v <= FY)
    print('duality: pairs', npairs, 'violations', dual)
    vals = [int(l.split()[1]) for l in open('../data/b212438.txt') if l.strip() and not l.startswith('#')]
    it = iter(vals)
    ok = bad = 0
    for n in range(4, 18):
        for k in range(4, 2 * n - 3):
            a = next(it)
            g = u(k, n)
            if g is None:
                continue
            ok += g == a
            bad += g != a
    print('A212438 cells ok', ok, 'bad', bad)
    for name, edge_val in (('b002840.txt', lambda v, f: u(v, f)),
                           ('b006445.txt', lambda v, f: None if u(v, f) is None else 2 * u(v, f) - S.get((v, f), S.get((f, v), 0)))):
        ref = bfile(name)
        ok = bad = skip = 0
        for e, want in ref.items():
            tot = 0
            for v in range(1, e + 2):
                f = e + 2 - v
                if v < 4 or f < 4:
                    continue
                g = edge_val(v, f)
                if g is None:
                    tot = None
                    break
                tot += g
            if tot is None:
                skip += 1
            elif tot == want:
                ok += 1
            else:
                bad += 1
                print(name, 'edges', e, 'got', tot, 'want', want)
        print(name, 'ok', ok, 'bad', bad, 'not covered', skip)
    known = bfile('b000944.txt')
    for n in range(4, VX + 1):
        if 2 * n - 4 > FY:
            continue
        a = sum(U.get((n, f), 0) for f in range(4, 2 * n - 3))
        tag = ('ok' if known[n] == a else f'MISMATCH (OEIS {known[n]})') if n in known else 'NEW'
        print(f'A000944({n}) = {a}  {tag}')


if __name__ == '__main__':
    main()
