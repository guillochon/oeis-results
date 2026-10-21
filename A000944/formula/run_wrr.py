"""All rotatory reflections W_rr mod one prime on the box v <= VX, f <= FY (M variables); cached .npy.
    python run_wrr.py VX FY prime_index
antipodal + even orders >= 4 (projective-plane quotients) + C_nh (odd n, symmetric mirror halves)."""
import os
import sys

import numpy as np

from disk import wrr_h
from ms import Ctx, no_throttle
from pipeline import halfturn, networks
from rotorefl import antipodal, even_order
from rp2 import maps_sphere_rp2
from t3 import PRIMES24, T3


def wrr_all(VX, FY, p):
    N = VX + FY
    out = np.zeros((VX + 1, FY + 1), dtype=np.int64)
    # projective-plane quotients: v_Q <= VX/2 + 1, f_Q <= FY/2 + 1, root face degree <= 2 e_Q
    T = T3(VX // 2 + 2, FY // 2 + 2, VX + FY + 6, p)
    P, Q = maps_sphere_rp2(T)
    h = Ctx(N // 2 + 1, p)
    KA = antipodal(h, P, Q, T)
    for a, b, c in KA.items():
        if 2 * a <= VX and 2 * b <= FY:
            out[2 * a, 2 * b] = (out[2 * a, 2 * b] + c) % p
    ctx = Ctx(N, p)
    nets = networks(ctx)
    _, Drev = halfturn(ctx, nets, want_rev=True)
    for n in range(4, N + 1, 2):
        KC = even_order(ctx, n, nets, Drev, P, Q, T)
        for t, ser in KC.items():
            for v, f, c in ser.items():
                if v <= VX and f <= FY:
                    out[v, f] = (out[v, f] + c) % p
        print('order', n, flush=True)
    out = (out + wrr_h(VX, FY, p)) % p
    return out


if __name__ == '__main__':
    no_throttle()
    VX, FY, i = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    p = PRIMES24[i]
    fn = f'wrr_{VX}_{FY}_{p}.npy'
    if not os.path.exists(fn):
        np.save(fn, wrr_all(VX, FY, p))
    print('done', VX, FY, p)
