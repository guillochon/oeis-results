import sys, glob, os, numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); from rules import codes, d4
P = {os.path.basename(f)[:-4]: np.load(f) for f in sorted(glob.glob('panels/*.npy'))}
def orb(w):
    B = np.array([(w >> (8 - q)) & 1 for q in range(9)]).reshape(3, 3)
    return min(int(''.join(map(str, M.ravel())), 2) for M in d4(B))
O = [orb(w) for w in range(512)]
cnt = Counter(); pan = {}
for k, A in P.items():
    c = Counter(O[w] for w in codes(A, 3, 3).ravel())   # one orientation is enough: panels are symmetric
    cnt.update(c)
    for o in c: pan.setdefault(o, set()).add(k)
tot = sum(cnt.values())
dens = np.mean([A.mean() for A in P.values()])
def expect(o):
    w = o; ones = bin(w).count('1'); size = len({x for x in range(512) if O[x] == o})
    return tot * size * dens ** ones * (1 - dens) ** (9 - ones)
rows = sorted(set(O), key=lambda o: cnt[o] / expect(o))
print('total windows', tot, 'density', round(dens, 3))
for o in rows[:16]:
    pic = ' '.join(''.join('#' if (o >> (8 - (3 * i + j))) & 1 else '.' for j in range(3)) for i in range(3))
    print(f'{pic}  used {cnt[o]:4d}  random-expect {expect(o):7.0f}  ratio {cnt[o]/expect(o):.3f}  in {len(pan.get(o,()))} panels')
