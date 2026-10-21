"""Sensed polyhedra S(v,f) (up to orientation-preserving isomorphism) = A239893, by Burnside:
S = 2 * (sum over rooted / (4e) + sum of rotation weights), computed mod several primes + CRT.

Usage: python sensed.py N [nprimes [VX FY]]   -> writes sensed_N.json {"v,f": S}
(with VX FY: only the cells v <= VX, f <= FY, which is all that combine_full.py VX FY reads)
"""
import json
import sys

from ms import Ctx, PRIMES, crt
from pipeline import halfturn, networks, rooted_cnets, rot_big


def rotation_total(ctx, nets):
    """Series in M-variables: sum over all rotation types of the weights (identity excluded)."""
    N = ctx.N
    tot = ctx.zero()
    K = halfturn(ctx, nets)
    for t in K:
        tot = tot + K[t]
    for k in range(3, N + 1):
        for tau, (pv, pf) in (('VF', (1, 1)), ('VV', (2, 0)), ('FF', (0, 2))):
            c = rot_big(ctx, k, tau, nets)
            for a, b, val in c.items():
                v, f = k * a + pv, k * b + pf
                if v + f <= N:
                    tot.a[v, f] = (tot.a[v, f] + val) % ctx.p
    return tot


def identity_weight(ctx, nets):
    """sum 1/|Aut| = rooted / (4e), in M-variables."""
    N = ctx.N
    m = rooted_cnets(ctx, nets)
    out = ctx.zero()
    for a, b, val in m.items():
        v, f = a + 2, b + 2
        if v + f <= N:
            e = v + f - 2
            out.a[v, f] = val * pow(4 * e, ctx.p - 2, ctx.p) % ctx.p
    return out


def sensed_mod(N, p):
    ctx = Ctx(N, p)
    nets = networks(ctx)
    W = identity_weight(ctx, nets) + rotation_total(ctx, nets)
    return (W * 2).a


def main():
    from ms import no_throttle
    no_throttle()
    N = int(sys.argv[1])
    npr = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    VX, FY = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (N, N)
    res = [sensed_mod(N, p) for p in PRIMES[:npr]]
    out = {}
    for v in range(4, N + 1):
        for f in range(4, min(N - v, FY) + 1):
            if v > VX:
                break
            val = crt([int(r[v, f]) for r in res], PRIMES[:npr])
            val2 = crt([int(r[v, f]) for r in res[:-1]], PRIMES[:npr - 1])
            assert val == val2, (v, f, 'CRT unstable: add primes')
            if val:
                out[f'{v},{f}'] = val
    json.dump(out, open(f'sensed_{N}.json', 'w'), indent=0)
    print('wrote', len(out), 'cells')


if __name__ == '__main__':
    main()
