"""Open pockets away from the center: sizes in the 19 panels, the panels holding the small ones, and how often
20 density-matched no-pinhole random designs per panel (null/, tools/null_gen.sh) have no off-center pocket under 5 or 9.
Run from attacks/vinciolo."""
import sys, os; sys.path.insert(0, "tools")
from growth_load import load
from collections import Counter
from scipy.ndimage import label
import numpy as np
P = load()
off = Counter(); cen = []; small = []
for name, A in P.items():
    n = A.shape[0]
    L, k = label(A == 0, structure=np.ones((3, 3)))
    c = sorted({(n - 1) // 2, n // 2})
    cl = {L[i, j] for i in c for j in c} - {0}
    for r in range(1, k + 1):
        s = int((L == r).sum())
        if r in cl: cen.append((name, s))
        else:
            off[s] += 1
            if s <= 8: small.append((name, s))
print('center regions:', cen)
print('smallest off-center sizes:', sorted(off.items())[:12])
print('off-center regions total:', sum(off.values()))
print('panels with each small off-center size:', sorted(small, key=lambda x: x[1]))
def offmin(A):
    n = A.shape[0]; L, k = label(A == 0, structure=np.ones((3, 3)))
    c = sorted({(n - 1) // 2, n // 2}); cl = {L[i, j] for i in c for j in c} - {0}
    s = [int((L == r).sum()) for r in range(1, k + 1) if r not in cl]
    return min(s) if s else 10**9
for T in (5, 9):
    joint = 1.0; rows = []
    for name, A in P.items():
        n = A.shape[0]; ok = tot = 0
        for l in open(os.path.join('null', name + '.txt')):
            B = np.pad(np.array([int(ch) for ch in l.strip()]).reshape(n - 2, n - 2), 1, constant_values=1)
            ok += offmin(B) >= T; tot += 1
        joint *= max(ok, 0.5) / tot; rows.append(ok)
    print(f'off-center min >= {T}: random designs passing per panel (of 20):', rows, ' joint ~', f'{joint:.1e}')
