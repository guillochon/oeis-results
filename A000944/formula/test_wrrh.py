import glob, json
from fractions import Fraction
from disk import wrr_h
from ms import frac_mod
from t3 import PRIMES24
p = PRIMES24[0]
VX, FY = 11, 18
W = wrr_h(VX, FY, p, verbose=True)
stats = {}
for fn in glob.glob('../brute/stats2_*.json'):
    for k, d in json.load(open(fn)).items():
        stats.setdefault(k, {}).update(d)
ok = bad = 0
for k, d in stats.items():
    v, f = map(int, k.split(','))
    if v > VX or f > FY:
        continue
    want = sum((Fraction(w) for t, w in d.items() if json.loads(t)[0] == 'rr' and json.loads(t)[-1] == 'h'), Fraction(0))
    if int(W[v, f]) == frac_mod(want, p): ok += 1
    else:
        bad += 1; print('mismatch', v, f, 'brute', want)
print('C_nh: ok', ok, 'bad', bad)
