"""Independent check of strategy (b): quintuples = 5-cliques of the tangency graph.

Spheres are integer row vectors u (sphere = u * C).  The Lorentz Gram matrix of
the starting Descartes configuration is G0 = 2I - J (self 1, tangent -1), and
every R_i preserves it, so spheres u, v are tangent iff u G0 v^T = -1.
Every set of 5 mutually tangent spheres found among the spheres of generation
<= K is counted, by the largest generation among its members.
"""
import sys
from collections import Counter
from itertools import combinations
from engine import strategy_b, reflect, IDENT, N

G0 = [[2 * (i == j) - 1 for j in range(N)] for i in range(N)]


def form(u, v):
    return sum(u[i] * G0[i][j] * v[j] for i in range(N) for j in range(N))


# R_i G0 R_i^T = G0
for i in range(N):
    R = reflect(IDENT, i)
    assert all(form(R[a], R[b]) == G0[a][b] for a in range(N) for b in range(N))

K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
quints, _, gen = strategy_b(K)
sph = [s for s, g in gen.items() if g <= K]
assert all(form(s, s) == 1 for s in sph)
idx = {s: k for k, s in enumerate(sph)}
import numpy as np
V = np.array(sph, dtype=np.int64)
VG = V @ np.array(G0, dtype=np.int64)
nbr = [set() for _ in sph]
for lo in range(0, len(sph), 2000):
    F = VG[lo:lo + 2000] @ V.T
    for a, b in zip(*np.nonzero(F == -1)):
        nbr[lo + a].add(int(b))

cnt = Counter()


def extend(clique, cand):
    if len(clique) == 5:
        cnt[max(gen[sph[c]] for c in clique)] += 1
        return
    for c in sorted(cand):
        if c > clique[-1]:
            extend(clique + [c], cand & nbr[c])


for a in range(len(sph)):
    extend([a], nbr[a])
print("cliques by gen", [cnt[g] for g in range(K + 1)])
print("BFS quintuples", quints[:K + 1])
print("max clique size check:", max(len(n) for n in nbr), "neighbors max")
