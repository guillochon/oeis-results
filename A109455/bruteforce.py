"""A109455: threshold functions of n variables up to permutations of the variables.

Method:
 1. Generate all *positive* threshold functions with sorted integer weights
    w1 >= ... >= wn >= 0, max weight <= B, every threshold t.
 2. Close under all n! variable permutations  -> labeled positive threshold functions
    (must equal A002078(n): 2,3,6,20,150,3287,244158).
 3. Close under all 2^n input negations        -> all threshold functions
    (must equal A000609(n): 2,4,14,104,1882,94572,15028134).
    Every generated function is a threshold function, so matching the known total
    proves the weight bound B was large enough (nothing missed).
 4. Count S_n-orbits two ways:
    (a) Burnside: (1/n!) * sum over cycle types of class_size * Fix(sigma)
    (b) direct: canonical form = min truth table over all n! permutations.

Truth table T of f: bit x of T is f(x), where bit i of x is variable x_i.
Usage: python A109455/bruteforce.py [max_n] [B]
"""
import sys
from itertools import combinations_with_replacement, permutations
from math import factorial
from collections import Counter

import numpy as np

A000609 = [2, 4, 14, 104, 1882, 94572, 15028134]
A002078 = [2, 3, 6, 20, 150, 3287, 244158]
A109455_LISTED = [2, 4, 10, 34, 178, 1720, 590440]


def point_perm(n, p):
    """Map on points x induced by sending variable i to position p[i]."""
    return [sum(((x >> i) & 1) << p[i] for i in range(n)) for x in range(1 << n)]


def apply_point_map(T, pm):
    """g(x) = f(pm[x]) for a numpy uint64 array of truth tables."""
    G = np.zeros_like(T)
    one = np.uint64(1)
    for x, px in enumerate(pm):
        G |= ((T >> np.uint64(px)) & one) << np.uint64(x)
    return G


def sorted_positive(n, B):
    pts = np.arange(1 << n)
    bits = ((pts[:, None] >> np.arange(n)) & 1)          # (2^n, n)
    out = set()
    for w in combinations_with_replacement(range(B, -1, -1), n):  # non-increasing weights
        s = bits @ np.array(w)                            # weighted sum at each point
        for t in range(0, int(s.max()) + 2):
            out.add(int(sum(1 << int(x) for x in np.nonzero(s >= t)[0])))
    return out


def all_threshold(n, B):
    base = np.array(sorted(sorted_positive(n, B)), dtype=np.uint64)
    perms = list(permutations(range(n)))
    pos = np.unique(np.concatenate([apply_point_map(base, point_perm(n, p)) for p in perms]))
    alls = [pos]
    for i in range(n):  # negate variable i: g(x) = f(x xor 2^i); close iteratively
        cur = np.concatenate(alls)
        alls = [np.unique(np.concatenate([cur, apply_point_map(cur, [x ^ (1 << i) for x in range(1 << n)])]))]
    return pos, alls[0]


def cycle_type_reps(n):
    """Yield (cycle_type, class_size, representative permutation as list p)."""
    counts = Counter(tuple(sorted(cycle_lengths(p), reverse=True)) for p in permutations(range(n)))
    for ct, size in sorted(counts.items()):
        p, start = [0] * n, 0
        for L in ct:
            for k in range(L):
                p[start + k] = start + (k + 1) % L
            start += L
        yield ct, size, p


def cycle_lengths(p):
    seen, out = set(), []
    for i in range(len(p)):
        if i not in seen:
            L, j = 0, i
            while j not in seen:
                seen.add(j); j = p[j]; L += 1
            out.append(L)
    return out


def burnside(n, F):
    total, rows = 0, []
    for ct, size, p in cycle_type_reps(n):
        fix = int(np.count_nonzero(apply_point_map(F, point_perm(n, p)) == F))
        rows.append((ct, size, fix))
        total += size * fix
    q, r = divmod(total, factorial(n))
    return q, r, rows


def direct_orbits(n, F):
    canon = F.copy()
    for p in permutations(range(n)):
        canon = np.minimum(canon, apply_point_map(F, point_perm(n, p)))
    return len(np.unique(canon))


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    for n in range(0, max_n + 1):
        pos, F = all_threshold(n, B)
        ok_pos = len(pos) == A002078[n]
        ok_all = len(F) == A000609[n]
        q, r, rows = burnside(n, F)
        line = (f"n={n}: positive={len(pos)} ({'OK' if ok_pos else 'MISMATCH'} A002078)  "
                f"all={len(F)} ({'OK' if ok_all else 'MISMATCH'} A000609)  "
                f"Burnside orbits={q} rem={r}  listed A109455={A109455_LISTED[n]}")
        if n <= 5 or "--direct" in sys.argv:
            d = direct_orbits(n, F)
            line += f"  direct orbits={d}"
        print(line, flush=True)
        if n == max_n:
            print("  Fix(sigma) by cycle type:")
            for ct, size, fix in rows:
                print(f"    {ct}: class size {size}, Fix {fix}")


if __name__ == "__main__":
    main()
