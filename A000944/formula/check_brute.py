"""Check the modular pipeline against the brute-force symmetry statistics (brute/stats_*.json)."""
import glob
import json
import sys
from fractions import Fraction

from ms import Ctx, PRIMES, frac_mod
from pipeline import halfturn, networks, rooted_cnets, rot_big


def load_stats():
    stats = {}
    for fn in glob.glob('../brute/stats_*.json'):
        for key, v in json.load(open(fn)).items():
            stats.setdefault(key, {}).update({t: Fraction(w) for t, w in v.items()})
    return stats


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    p = PRIMES[0]
    ctx = Ctx(N, p)
    stats = load_stats()
    nets = networks(ctx)
    bad = checked = 0
    m = rooted_cnets(ctx, nets)
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v + f - 4 <= N:
            e = v + f - 2
            checked += 1
            if m[v - 2, f - 2] != frac_mod(4 * e * d['["id"]'], p):
                bad += 1
                print('rooted mismatch', v, f)
    for k in range(3, N):
        for tau, (pv, pf) in (('VF', (1, 1)), ('VV', (2, 0)), ('FF', (0, 2))):
            c = rot_big(ctx, k, tau, nets)
            pol = ''.join(sorted(tau))
            for key, d in stats.items():
                v, f = map(int, key.split(','))
                want = frac_mod(d.get(json.dumps(['rot', k, pol]), 0), p)
                got = 0
                if (v - pv) % k == 0 and (f - pf) % k == 0:
                    a, b = (v - pv) // k, (f - pf) // k
                    if a + b > N:
                        continue
                    got = c[a, b]
                checked += 1
                if got != want:
                    bad += 1
                    print('rot mismatch', k, tau, v, f)
    K = halfturn(ctx, nets)
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v + f > N:
            continue
        for t in K:
            pol = ''.join(sorted(t))
            want = frac_mod(d.get(json.dumps(['rot', 2, pol]), 0), p)
            checked += 1
            if K[t][v, f] != want:
                bad += 1
                print('half-turn mismatch', t, v, f, 'want', d.get(json.dumps(['rot', 2, pol]), 0))
    print('checked', checked, 'bad', bad)


if __name__ == '__main__':
    main()
