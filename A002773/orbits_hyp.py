"""Second, independent count of loopless non-simple matroids on n <= 10 elements.

Every loopless matroid N on n elements is a parallel extension of its simplification M
(simple, k points): replace point i by a parallel class of size w_i >= 1, sum w = n.
N ~ N' iff (M, w) ~ (M', w'), i.e. M = M' and w' = w o g for some g in Aut(M).
So #loopless(n) = sum over simple M on k <= n points of #(Aut(M)-orbits on compositions of n).

Method (deliberately different from experiments/deep_dives_nov2/deep_alg/orb.c, which used a
rank-function backtrack and Burnside's lemma):
  * Aut(M) = permutations preserving the set of hyperplanes (a simple matroid is determined by
    its hyperplanes).
  * A strong generating set is found level by level (one automorphism per coset of the point
    stabilizer chain), pruned by pairwise hyperplane-size profiles; leaves are checked exactly.
  * Orbits on compositions are counted explicitly by union-find under the generators.

Also: labeled counts of all matroids on k <= 9 points by rank, as sum of k!/|Aut(M)|.

Usage: python orbits_hyp.py <matroids09_rankLine> [workers]
"""
import os
import sys
from itertools import combinations
from math import factorial
from multiprocessing import Pool

NMAX = 10


def compositions(s, k):
    # all k-tuples of positive integers summing to s
    out = []
    for cuts in combinations(range(1, s), k - 1):
        b = (0,) + cuts + (s,)
        out.append(tuple(b[i + 1] - b[i] for i in range(k)))
    return out


COMPS = {(s, k): compositions(s, k) for k in range(1, 10) for s in range(k, NMAX + 1)}


def hyperplanes(line, k, r):
    full = (1 << k) - 1
    hs = []
    for X in range(1 << k):
        if ord(line[X]) - 48 != r - 1:
            continue
        if all(ord(line[X | 1 << e]) - 48 == r for e in range(k) if not X >> e & 1):
            hs.append(X)
    return hs


def apply(p, X):
    Y = 0
    while X:
        b = X & -X
        Y |= 1 << p[b.bit_length() - 1]
        X ^= b
    return Y


def find_aut(k, hyps, hypset, inv2, prefix):
    """One automorphism extending `prefix` (images of 0..len-1, in order), or None.
    Pruning: inv2[a][b] = sorted sizes of the hyperplanes containing a and b must be preserved."""
    p = [None] * k
    used = [False] * k

    def fits(m, t):
        return all(inv2[j][m] == inv2[p[j]][t] for j in range(m)) and inv2[m][m] == inv2[t][t]

    for m, t in prefix:
        if used[t] or not fits(m, t):
            return None
        p[m] = t
        used[t] = True

    def bt(m):
        if m == k:
            return all(apply(p, H) in hypset for H in hyps)
        for t in range(k):
            if used[t] or not fits(m, t):
                continue
            p[m] = t
            used[t] = True
            if bt(m + 1):
                return True
            p[m] = None
            used[t] = False
        return False

    return list(p) if bt(len(prefix)) else None


def generators(k, hyps):
    hypset = set(hyps)
    sizes = {H: bin(H).count("1") for H in hyps}
    inv2 = [[tuple(sorted(sizes[H] for H in hyps if H >> a & 1 and H >> b & 1)) for b in range(k)]
            for a in range(k)]
    gens = []
    order = 1  # |Aut| = product of the basic orbit lengths |i^(G_{0..i-1})|
    for i in range(k):
        base = [(j, j) for j in range(i)]
        orbit = 1
        for t in range(i + 1, k):
            g = find_aut(k, hyps, hypset, inv2, base + [(i, t)])
            if g is not None:
                gens.append(g)
                orbit += 1
        order *= orbit
    return gens, order


def orbit_count(comps, gens, k):
    idx = {c: n for n, c in enumerate(comps)}
    par = list(range(len(comps)))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for g in gens:
        for c, n in idx.items():
            d = [0] * k
            for i in range(k):
                d[g[i]] = c[i]
            a, b = f(n), f(idx[tuple(d)])
            if a != b:
                par[a] = b
    return sum(1 for n in range(len(comps)) if f(n) == n)


def work(chunk):
    tot = {}
    nsimple = {}
    labeled = {}
    for line in chunk:
        if line == "NULL":
            continue  # empty matroid: contributes only to s = 0
        k = len(line).bit_length() - 1
        r = ord(line[-1]) - 48
        # hyperplanes determine any matroid, so this also gives Aut(M) for non-simple M
        hyps = hyperplanes(line, k, r)
        gens, order = generators(k, hyps)
        assert factorial(k) % order == 0
        labeled[k, r] = labeled.get((k, r), 0) + factorial(k) // order
        if any(line[1 << e] != "1" for e in range(k)):
            continue
        if any(line[(1 << a) | (1 << b)] != "2" for a in range(k) for b in range(a + 1, k)):
            continue
        nsimple[k, r] = nsimple.get((k, r), 0) + 1
        for s in range(k, NMAX + 1):
            o = orbit_count(COMPS[s, k], gens, k)
            tot[s, r] = tot.get((s, r), 0) + o
    return tot, nsimple, labeled


def main():
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    path = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    lines = open(path).read().split()
    chunks = [lines[i::200] for i in range(200)]
    tot, nsimple, labeled = {}, {}, {(0, 0): 1}
    with Pool(workers) as pool:
        for t, ns, lab in pool.imap_unordered(work, chunks):
            for key, v in lab.items():
                labeled[key] = labeled.get(key, 0) + v
            for key, v in t.items():
                tot[key] = tot.get(key, 0) + v
            for key, v in ns.items():
                nsimple[key] = nsimple.get(key, 0) + v
    tot[0, 0] = 1
    for k in range(10):
        print("simple k=%d:" % k, *[nsimple.get((k, r), int(k == r == 0)) for r in range(10)])
    for s in range(NMAX + 1):
        row = [tot.get((s, r), 0) for r in range(NMAX + 1)]
        print("loopless-from-k<=9 s=%d:" % s, *row, " sum", sum(row))
    for k in range(10):
        row = [labeled.get((k, r), 0) for r in range(k + 1)]
        print("labeled all n=%d:" % k, *row, " sum", sum(row))


if __name__ == "__main__":
    main()
