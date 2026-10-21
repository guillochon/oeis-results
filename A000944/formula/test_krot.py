"""sum_d Kvertex_n(v,f,d) must equal 2 c_VV + c_VF of rotations of order n (rot_big)."""
import sys
import numpy as np
from disk import kface_rotational
from ms import Ctx
from pipeline import networks, rot_big
from t3 import PRIMES24
p = PRIMES24[0]
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 18
for n in (3, 5):
    Kf = kface_rotational(n, NC, p)
    tot = np.zeros((NC + 1, NC + 1), dtype=np.int64)
    for d, arr in Kf.items():
        tot = (tot + arr.T) % p          # back to pole-vertex form, indexed [v, f]
    ctx = Ctx(2 * NC, p)
    nets = networks(ctx)
    want = np.zeros_like(tot)
    for tau, (pv, pf), mult in (('VV', (2, 0), 2), ('VF', (1, 1), 1)):
        c = rot_big(ctx, n, tau, nets)
        for a, b, val in c.items():
            v, f = n * a + pv, n * b + pf
            if v <= NC and f <= NC:
                want[v, f] = (want[v, f] + mult * val) % p
    ok = int(np.sum((tot == want)))
    bad = [(v, f) for v in range(NC + 1) for f in range(NC + 1) if tot[v, f] != want[v, f]]
    print('n', n, 'cells equal', ok, 'mismatches', bad[:10])
