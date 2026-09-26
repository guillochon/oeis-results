"""Check the conjectured closed forms for A236553 and A227867 against a direct count for every modulus
n <= MAX and against the OEIS data.

A236553(n): solutions mod n of  a^2 - b^2 + c^2 + d^2 = 1,  2ab = 2ac = 2ad = 0   (i^2 = j^2 = +1)
A227867(n): solutions mod n of  a^2 - b^2 - c^2 - d^2 = 1,  2ab = 2ac = 2ad = 0   (i^2 = j^2 = -1)
Usage: python verify.py [MAX]
"""
import sys

import numpy as np
from sympy import factorint


def count(n, signs):
    """signs = (sb, sc, sd): coefficients of b^2, c^2, d^2 in the first equation."""
    xs = np.arange(n, dtype=np.int64)
    sq = (xs * xs) % n
    cache, total = {}, 0
    for a in range(n):
        allowed = (2 * a * xs) % n == 0
        key = int(allowed.sum())  # the allowed set is a subgroup, fixed by its size
        if key not in cache:
            acc = np.zeros(n, dtype=np.int64)
            acc[0] = 1
            for s in signs:  # convolve in the distribution of s*x^2 over allowed x
                h = np.bincount((s * sq[allowed]) % n, minlength=n)
                new = np.zeros(n, dtype=np.int64)
                for r in np.nonzero(h)[0]:
                    new += h[r] * np.roll(acc, r)
                acc = new
            cache[key] = acc
        total += int(cache[key][(1 - a * a) % n])
    return total


def local_A236553(p, k):
    if p == 2:
        return {1: 8, 2: 64}.get(k, 2 ** (2 * k + 2) + 32)
    return p ** (2 * k - 1) * (p + 1) + 2


def local_A227867(p, k):
    if p == 2:
        return {1: 8}.get(k, 32)
    return p ** (2 * k - 1) * (p + 1) + 2


def closed(n, local):
    out = 1
    for p, k in factorint(n).items():
        out *= local(p, k)
    return out


OEIS = {
    "A236553": [1, 8, 14, 64, 32, 112, 58, 288, 110, 256, 134, 896, 184, 464, 448, 1056, 308, 880, 382, 2048, 812,
                1072, 554, 4032, 752, 1472, 974, 3712, 872, 3584, 994, 4128, 1876, 2464, 1856, 7040, 1408, 3056,
                2576, 9216, 1724, 6496, 1894, 8576, 3520, 4432, 2258, 14784, 2746],
    "A227867": [1, 8, 14, 32, 32, 112, 58, 32, 110, 256, 134, 448, 184, 464, 448, 32, 308, 880, 382, 1024, 812, 1072,
                554, 448, 752, 1472, 974, 1856, 872, 3584, 994, 32, 1876, 2464, 1856, 3520, 1408, 3056, 2576, 1024,
                1724, 6496, 1894, 4288, 3520, 4432, 2258, 448, 2746, 6016, 4312, 5888],
}

if __name__ == "__main__":
    MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    for name, signs, local in [("A236553", (-1, 1, 1), local_A236553), ("A227867", (-1, -1, -1), local_A227867)]:
        bad = [n for n in range(1, MAX + 1) if count(n, signs) != closed(n, local)]
        data_bad = [n for n, v in enumerate(OEIS[name], 1) if v != closed(n, local)]
        print(f"{name}: direct count vs closed form, n=1..{MAX}: {'all match' if not bad else 'MISMATCH at ' + str(bad[:10])}")
        print(f"{name}: OEIS data ({len(OEIS[name])} terms) vs closed form: {'all match' if not data_bad else 'MISMATCH at ' + str(data_bad)}")
