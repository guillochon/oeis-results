"""Completeness check for strategy (b) with geometry only: every set of five mutually tangent objects
(distance between centers = r1 + r2; sphere-plate: |z| + r = 1) among the objects of generation <= K,
counted by its youngest member, must equal the quintuple counts of the search.

Uses the sphere positions and generations of verify_geom.strategy_b (no integer vectors anywhere).
"""
import sys
from collections import Counter

import numpy as np

from verify_geom import strategy_b, Sph

K = int(sys.argv[1]) if len(sys.argv) > 1 else 7
_, quints, st = strategy_b(K, check=False)
idx = [i for i, g in enumerate(st.gen) if g <= K]
gen = np.array([st.gen[i] for i in idx])
sph = [i for i in idx if isinstance(st.obj[i], Sph)]
pl = [i for i in idx if not isinstance(st.obj[i], Sph)]
pos = {i: k for k, i in enumerate(idx)}
C = np.array([st.obj[i].c for i in sph]); R = np.array([st.obj[i].r for i in sph])
nbr = [set() for _ in idx]
B = 2000
for lo in range(0, len(sph), B):
    d = np.linalg.norm(C[lo:lo + B, None, :] - C[None, :, :], axis=2)
    s = R[lo:lo + B, None] + R[None, :]
    a, b = np.nonzero(np.abs(d - s) < 1e-9 * np.maximum(1, s))
    for x, y in zip(a, b):
        if lo + x != y:
            nbr[pos[sph[lo + x]]].add(pos[sph[y]])
for p in pl:
    P = st.obj[p]
    t = np.abs(np.abs(C @ P.n - P.h) - R) < 1e-9 * np.maximum(1, R)
    for k in np.nonzero(t)[0]:
        nbr[pos[p]].add(pos[sph[k]]); nbr[pos[sph[k]]].add(pos[p])
# the two plates touch each other (at infinity) and are part of the starting quintuple
nbr[pos[pl[0]]].add(pos[pl[1]]); nbr[pos[pl[1]]].add(pos[pl[0]])

cnt = Counter()


def extend(cl, cand):
    if len(cl) == 5:
        cnt[int(max(gen[c] for c in cl))] += 1
        return
    for c in cand:
        if c > cl[-1]:
            extend(cl + [c], cand & nbr[c])


for a in range(len(idx)):
    extend([a], nbr[a])
print("geometric 5-cliques by youngest member:", [cnt[g] for g in range(K + 1)])
print("quintuples from the search:            ", quints[:K])
