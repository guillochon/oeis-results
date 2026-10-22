"""Sensed counts S(v,f) mod one prime, cached, so the primes can run in parallel; then combine.

    python sensed_p.py N i VX FY     -> sensed_N_p.npy (S mod PRIMES[i] on v <= VX, f <= FY, v + f <= N)
    python sensed_p.py N combine nprimes VX FY   -> sensed_N.json (CRT, checked with one prime fewer)
"""
import json
import os
import sys

import numpy as np

from ms import PRIMES, crt, no_throttle
from sensed import sensed_mod

if __name__ == '__main__':
    no_throttle()
    N = int(sys.argv[1])
    if sys.argv[2] != 'combine':
        i, VX, FY = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        p = PRIMES[i]
        fn = f'sensed_{N}_{p}.npy'
        if not os.path.exists(fn):
            a = sensed_mod(N, p)
            np.save(fn, a[:VX + 1, :FY + 1])
        print('done', N, p)
    else:
        npr, VX, FY = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        res = [np.load(f'sensed_{N}_{p}.npy') for p in PRIMES[:npr]]
        out = {}
        for v in range(4, VX + 1):
            for f in range(4, min(N - v, FY) + 1):
                val = crt([int(r[v, f]) for r in res], PRIMES[:npr])
                assert val == crt([int(r[v, f]) for r in res[:-1]], PRIMES[:npr - 1]), (v, f, 'CRT unstable')
                if val:
                    out[f'{v},{f}'] = val
        json.dump(out, open(f'sensed_{N}.json', 'w'), indent=0)
        print('wrote', len(out), 'cells')
