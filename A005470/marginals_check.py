"""Marginals of the formula cycle index Z_K (kf_N_q.json rows) against exhaustive-generation data.

Summing the rows' weights over all cycle types counts polyhedra up to isomorphism (Burnside), so mod
each prime q the rows must reproduce:
  A000944  by vertices (plantri, n <= 18)        A002840  by edges (e <= 35)
  A049337  by vertices and edges (n <= 14)       A000109  triangulations, e = 3n - 6 (n <= 24)
Usage: python marginals_check.py N NPRIMES
"""
import json
import os
import sys
from collections import defaultdict

import tower_p as T

A000944 = [0, 0, 0, 1, 2, 7, 34, 257, 2606, 32300, 440564, 6384634, 96262938, 1496225352, 23833988129,
           387591510244, 6415851530241, 107854282197058]                                   # n = 1..
A002840 = [1, 0, 1, 2, 2, 4, 12, 22, 58, 158, 448, 1342, 4199, 13384, 43708, 144810, 485704, 1645576,
           5623571, 19358410, 67078828, 233800162, 819267086, 2884908430, 10204782956, 36249143676,
           129267865144, 462669746182, 1661652306539, 5986979643542]                       # e = 6..
A000109 = [1, 1, 1, 2, 5, 14, 50, 233, 1249, 7595, 49566, 339722, 2406841, 17490241, 129664753,
           977526957, 7475907149, 57896349553, 453382272049, 3585853662949, 28615703421545,
           230214772514389]                                                                 # n = 3..


def a049337():
    f = os.path.expanduser('~/oeis-mcp/data/oeisdata/files/A049/b049337.txt')
    vals = [int(l.split()[1]) for l in open(f) if l.strip() and not l.startswith('#')]
    tri, i, n = {}, 4, 4        # the b-file's rows n = 1..3 take 4 entries in all; row n has C(n,2)+1
    while i + n * (n - 1) // 2 + 1 <= len(vals):
        for e in range(n * (n - 1) // 2 + 1):
            tri[(n, e)] = vals[i + e]
        i += n * (n - 1) // 2 + 1
        n += 1
    return tri


def main():
    N, npr = int(sys.argv[1]), int(sys.argv[2])
    tri = a049337()
    bad = 0
    for q in T.PRIMES24[:npr]:
        byn, bye, byne = defaultdict(int), defaultdict(int), defaultdict(int)
        for n, vp, eo, w in json.load(open(f'kf_{N}_{q}.json')):
            e = sum(k for k, _ in eo)
            byn[n] += w
            bye[e] += w
            byne[(n, e)] += w
        checks = []
        checks += [('A000944', n, byn[n], v) for n, v in enumerate(A000944, 1) if n >= 4 and n <= N]
        checks += [('A002840', e, bye[e], v) for e, v in enumerate(A002840, 6) if (2 * e) // 3 + 2 <= N]
        checks += [('A049337', k, byne[k], v) for k, v in tri.items() if k[0] >= 4 and k[0] <= N]
        checks += [('A000109', n, byne[(n, 3 * n - 6)], v) for n, v in enumerate(A000109, 3)
                   if 4 <= n <= N]
        cnt = defaultdict(lambda: [0, 0])
        for a, k, got, want in checks:
            cnt[a][1] += 1
            if got % q == want % q:
                cnt[a][0] += 1
            else:
                bad += 1
                print('  MISMATCH', a, k, q)
        print(q, '  '.join(f'{a} {o}/{t}' for a, (o, t) in cnt.items()), flush=True)
    print('ALL OK' if not bad else f'{bad} MISMATCHES')


if __name__ == '__main__':
    main()
