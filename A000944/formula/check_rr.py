"""Check rotatory reflections (antipodal and even order >= 4) against brute/stats2_*.json."""
import glob
import json
import sys
from fractions import Fraction

from ms import Ctx, PRIMES, frac_mod
from pipeline import halfturn, networks
from rotorefl import antipodal, even_order
from rp2 import maps_sphere_rp2
from t3 import T3

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
p = PRIMES[0]
stats = {}
for fn in glob.glob('../brute/stats2_*.json'):
    for key, d in json.load(open(fn)).items():
        stats.setdefault(key, {}).update(d)
vmax = max(int(k.split(',')[0]) for k in stats)
ctx = Ctx(N, p)
nets = networks(ctx)
K, Drev = halfturn(ctx, nets, want_rev=True)
T = T3(N + 2, N + 2, 2 * N + 6, p)
P, Q = maps_sphere_rp2(T)
ok = bad = 0
for n in range(4, N + 1, 2):
    KC = even_order(ctx, n, nets, Drev, P, Q, T)
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v + f > N:
            continue
        for t, ser in KC.items():
            want = Fraction(d.get(json.dumps(['rr', n, ''.join(sorted(t))]), '0'))
            if ser[v, f] == frac_mod(want, p):
                ok += 1
            else:
                bad += 1
                print('mismatch order', n, t, v, f, 'brute', want)
print('even orders: ok', ok, 'bad', bad)
