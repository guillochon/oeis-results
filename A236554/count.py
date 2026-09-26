"""Count involutions of the quaternion ring over Z/(2^n) with i^2 = j^2 = 1, straight from the definition.

Solutions mod N = 2^n of:  a^2 - b^2 + c^2 + d^2 = 1,  2ab = 2ac = 2ad = 0.
For each a, the allowed b, c, d form the set {x : 2ax = 0 mod N}; the number of triples hitting
the target 1 - a^2 is a cyclic convolution of square-value histograms (numpy FFT-free, exact ints).
Usage: python count.py [max_n]
"""
import sys

import numpy as np


def count(n):
    N = 1 << n
    xs = np.arange(N, dtype=np.int64)
    sq = (xs * xs) % N
    total = 0
    cache = {}
    for a in range(N):
        allowed = (2 * a * xs) % N == 0
        key = int(allowed.sum())  # allowed set depends only on gcd(2a, N)
        if key not in cache:
            h = np.bincount(sq[allowed], minlength=N)           # counts of x^2 for x allowed
            hneg = np.bincount((-sq[allowed]) % N, minlength=N)  # counts of -x^2
            cd = np.zeros(N, dtype=np.int64)                    # c^2 + d^2
            for r in np.nonzero(h)[0]:
                cd += h[r] * np.roll(h, r)
            bcd = np.zeros(N, dtype=np.int64)                   # c^2 + d^2 - b^2
            for r in np.nonzero(hneg)[0]:
                bcd += hneg[r] * np.roll(cd, r)
            cache[key] = bcd
        total += int(cache[key][(1 - a * a) % N])
    return total


if __name__ == "__main__":
    known = [8, 64, 288, 1056, 4128, 16416, 65568, 262176]
    for n in range(1, (int(sys.argv[1]) if len(sys.argv) > 1 else 14) + 1):
        c = count(n)
        formula = 2 ** (2 * n + 2) + 32 if n >= 3 else None
        tag = "OEIS ok" if n <= len(known) and known[n - 1] == c else ("OEIS MISMATCH" if n <= len(known) else "new")
        print(f"n={n:2d}  a(n)={c:<12d} 2^(2n+2)+32={formula}  {'match' if c == formula else ''}  [{tag}]")
