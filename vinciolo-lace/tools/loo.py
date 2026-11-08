import sys, glob, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__)); from rules import table, violations, codes, d4
P = {os.path.basename(f)[:-4]: np.load(f) for f in sorted(glob.glob('panels/*.npy'))}
for a, b in ((3, 3), (3, 4), (4, 4)):
    ok = table(P.values(), a, b)
    loo = {k: violations(A, a, b, table([B for j, B in P.items() if j != k], a, b)) for k, A in P.items()}
    print(f'{a}x{b}: all {len(P)} panels allow {ok.sum()} of {ok.size}; leave-one-out: '
          f'{sum(v == 0 for v in loo.values())} of {len(P)} pass, violations {sorted(loo.values())}')
ok = table(P.values(), 3, 3)
F = [w for w in range(512) if not ok[w]]
seen = set(); reps = []
for w in F:
    B = np.array([(w >> (8 - q)) & 1 for q in range(9)]).reshape(3, 3)
    o = min(int(''.join(map(str, M.ravel())), 2) for M in d4(B))
    if o not in seen: seen.add(o); reps.append(o)
print(len(F), 'forbidden 3x3 blocks in', len(reps), 'symmetry classes:')
for w in reps:
    print(' '.join(''.join('#' if (w >> (8 - (3 * i + j))) & 1 else '.' for j in range(3)) for i in range(3)))
