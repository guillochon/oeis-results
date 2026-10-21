import numpy as np
from disk import Grammar, rooted_by_anchors
from refl import wref
from t3 import PRIMES24, T3, cnets_by_rootface
p = PRIMES24[0]
# 1. refl mode reproduces refl.wref
VX, FY = 9, 14
G = Grammar('refl', VX, FY, p)
co = rooted_by_anchors(G)
W = G.R.zero()
for A in range(3, len(co)):
    W = (W + co[A] * pow(4 * A, p - 2, p)) % p
print('refl mode == refl.wref:', np.array_equal(W, wref(VX, FY, p)))
# 2. cone mode: [alpha^A] * x  == rooted c-nets with root vertex of degree A (dual of root face degree)
VX, FY = 12, 12
G = Grammar('cone', VX, FY, p)
co = rooted_by_anchors(G)
T = T3(FY + 1, VX + 1, VX + 1, p)
c, _ = cnets_by_rootface(T)
bad = ok = 0
for A in range(3, VX + 1):
    for v in range(4, VX + 1):
        for f in range(4, FY + 1):
            got = int(co[A][v - 1, f]) if v - 1 <= VX else None
            want = int(c[f, v, A]) if f <= T.X and v <= T.Y and A <= T.Z else None
            if got is None or want is None:
                continue
            if got == want: ok += 1
            else:
                bad += 1
                if bad < 10: print('cone mismatch A v f', A, v, f, got, want)
print('cone mode vs root-vertex degree: ok', ok, 'bad', bad)
