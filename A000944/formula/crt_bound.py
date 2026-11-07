"""Proved coefficient bound for the CRT reconstructions of S(v,f) and U(v,f).

combine_full.py and sensed_p.py recover integers from residues by the symmetric CRT lift, which returns
the true value only if the product P of the primes exceeds 2 max|value|. This script proves that
condition for every cell, with exact integer arithmetic throughout (it forces fastmul.EXACT).

  0 <= U(v,f) <= S(v,f) <= R(v,f) <= B(v,f)

  R = rooted c-nets with v vertices and f faces. Every polyhedron (sensed type) has at least one rooting,
      and different sensed types give different rooted maps, so S <= R; U <= S since each unsensed type
      is one or two sensed types.
  B = rooted nonseparable planar maps with v vertices and f faces (Brown-Tutte):
      B = (2i+j-2)! (2j+i-2)! / (i! j! (2i-1)! (2j-1)!) with i = v-1, j = f-1. A c-net is nonseparable,
      so R <= B.

Step 1: compute R mod enough 26-bit primes that their product exceeds B (R >= 0, so the least
nonnegative lift is R itself); R is then known exactly. Step 2: check 2 R(v,f) < P for the primes that
produced each stored table, and S <= R, U <= S on every cell.

Usage: python crt_bound.py [workers]     (about 2 minutes per prime in exact mode)
"""
import json
import math
import os
import sys

os.environ['FASTMUL_EXACT'] = '1'

from ms import PRIMES as SPRIMES  # noqa: E402
from t3 import PRIMES24  # noqa: E402

N = 146                      # v + f <= 146 covers the 50 x 96 box
U_FILE, U_PRIMES = 'Ufull_50_96.json', PRIMES24[:10]
S_FILE, S_PRIMES = 'sensed_146.json', SPRIMES[:10]


def brown_tutte(v, f):
    i, j = v - 1, f - 1
    F = math.factorial
    return F(2 * i + j - 2) * F(2 * j + i - 2) // (F(i) * F(j) * F(2 * i - 1) * F(2 * j - 1))


def is_prime(n):
    if n < 2:
        return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17):      # deterministic for n < 3.4e14
        x = pow(a, d, n)
        if x in (1, n - 1) or a % n == 0:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_below(top, count):
    out, n = [], top - 1
    while len(out) < count:
        if is_prime(n):
            out.append(n)
        n -= 1
    return out


def rooted_mod(p):
    from ms import Ctx, no_throttle
    from pipeline import networks, rooted_cnets
    no_throttle()
    ctx = Ctx(N, p)
    m = rooted_cnets(ctx, networks(ctx))
    return {(a + 2, b + 2): int(m[a, b]) for a in range(N + 1) for b in range(N + 1 - a)
            if a + 2 >= 4 and b + 2 >= 4 and a + b + 4 <= N}


def crt_nonneg(res, primes):
    x, m = 0, 1
    for r, p in zip(res, primes):
        x += m * ((r - x) * pow(m, -1, p) % p)
        m *= p
    return x


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    S = {tuple(map(int, k.split(','))): v for k, v in json.load(open(S_FILE)).items()}
    U = {tuple(map(int, k.split(','))): v for k, v in json.load(open(U_FILE)).items()}
    cells = sorted(set(S) | set(U))
    B = {c: brown_tutte(*c) for c in cells}
    need = 2 * max(B.values())
    primes, P = [], 1
    for p in primes_below(1 << 26, 40):
        if P > need:
            break
        primes.append(p)
        P *= p
    print(f'{len(cells)} cells; max B has {max(B.values()).bit_length()} bits; {len(primes)} primes, '
          f'{P.bit_length()} bits', flush=True)
    from multiprocessing import Pool
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    with Pool(workers) as pool:
        res = pool.map(rooted_mod, primes)
    R = {c: crt_nonneg([r.get(c, 0) for r in res], primes) for c in cells}
    bad = 0
    for c in cells:
        if not R[c] <= B[c]:
            bad += 1
            print('R > B at', c)
        if c in S and not S[c] <= R[c]:
            bad += 1
            print('S > R at', c)
        if c in U and not U[c] <= S.get(c, 0):
            bad += 1
            print('U > S at', c)
    PU, PS = math.prod(U_PRIMES), math.prod(S_PRIMES)
    rU = max(R[c] for c in U)
    rS = max(R[c] for c in S)
    print(f'U table: max R over its cells {rU.bit_length()} bits, CRT modulus {PU.bit_length()} bits, '
          f'2R < P: {2 * rU < PU}')
    print(f'S table: max R over its cells {rS.bit_length()} bits, CRT modulus {PS.bit_length()} bits, '
          f'2R < P: {2 * rS < PS}')
    ok = not bad and 2 * rU < PU and 2 * rS < PS
    print('CRT bound', 'PROVED for every cell' if ok else 'FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
