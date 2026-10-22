"""A000109 (simplicial polyhedra by vertices) = (A253882 + achiral) / 2, achiral = 2 (W_ref + W_rr) on the
triangulation diagonal f = 2n - 4.

A253882 (orientation-preserving count) has a closed formula (Howroyd); W_ref comes from triang.py (one
variable, fast); W_rr is read from run_wrr.py's arrays wrr_VX_FY_p.npy with FY = 2 VX - 4.

    python a109.py N nprimes     ->  b000109 terms for n = 4..N, checks, writes ../derived/a109_N.json
"""
import json
import sys
from math import comb, gcd

import numpy as np

from ms import crt
from t3 import PRIMES24
from triang import wref_triang


def phi(n):
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


def a253882(n):
    """Howroyd's formula (A253882, Mar 02 2021): triangulations up to orientation-preserving isomorphism."""
    from fractions import Fraction
    if n < 3:
        return 0
    t = Fraction(2 * comb(4 * (n - 3) + 1, n - 3), (n - 2) * (3 * n - 7))
    t += Fraction(3 * sum(phi(d) * comb(4 * ((n - 2) // d), (n - 2) // d) for d in range(2, n - 1) if (n - 2) % d == 0), 4)
    if n % 2 == 1:
        s = (n - 3) // 2
        t += Fraction(3 * comb(4 * s, s) * (2 * s + 1), 3 * s + 1)
    if n % 3 == 1:
        s = (n - 4) // 3
        t += Fraction(8 * comb(4 * s, s) * (4 * s + 1), 3 * s + 1)
    if n % 3 == 0:
        s = (n - 3) // 3
        t += 2 * comb(4 * s, s)
    t /= 6 * (n - 2)
    assert t.denominator == 1
    return int(t)


def main():
    N = int(sys.argv[1])
    npr = int(sys.argv[2])
    VX, FY = N, 2 * N - 4
    primes = PRIMES24[:npr]
    res = {n: [] for n in range(4, N + 1)}
    for p in primes:
        wr = wref_triang(N, p)
        wq = np.load(f'wrr_{VX}_{FY}_{p}.npy')
        for n in range(4, N + 1):
            res[n].append(2 * (int(wr[n]) + int(wq[n, 2 * n - 4])) % p)
    out = {}
    for n in range(4, N + 1):
        ach = crt(res[n], primes)
        assert ach == crt(res[n][:-1], primes[:-1]) and ach >= 0, ('CRT unstable / negative', n, ach)
        s = a253882(n)
        assert (s + ach) % 2 == 0 and ach <= s, ('parity / size', n)
        out[n] = (s + ach) // 2
    known = {int(l.split()[0]): int(l.split()[1]) for l in open('../derived/b000109.txt')}
    agree = [n for n in out if n in known and known[n] == out[n]]
    diff = [n for n in out if n in known and known[n] != out[n]]
    print(f'agree with b000109 (n <= 40, from the A000944 run): {len(agree)} terms; disagree: {diff}')
    json.dump({str(n): v for n, v in out.items()}, open(f'../derived/a109_{N}.json', 'w'), indent=0)
    for n in list(out)[-3:]:
        print(f'A000109({n}) = {out[n]}')


if __name__ == '__main__':
    main()
