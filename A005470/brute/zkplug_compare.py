"""Compare plantri_zk output (zkplug.c) with the brute-force zk*.json (exact) or the formula rows kf_N_q.json
(mod each prime), cell by cell.
Usage: python zkplug_compare.py n out_part1.txt [...] --exact zk.json
       python zkplug_compare.py n out_part1.txt [...] --kf N NPRIMES
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, '..')


def load_plug(n, files):
    tot = defaultdict(Fraction)
    for fn in files:
        for line in open(fn):
            s = line.replace('>>planar_code<<', '').split()
            if len(s) != 3:
                continue
            key, g, c = s[0], int(s[1]), int(s[2])
            if key == 'ID':
                vp = (1,) * n
                eo = None
            elif '|' in key:
                a, b = key.split('|')
                vp = tuple(int(x) for x in a.split(','))
                eo = tuple(tuple(int(y) for y in x.split('.')) for x in b.split(','))
            else:
                continue
            tot[(vp, eo)] += Fraction(c, g)
    # identity: edges (1,0)^e with e unknown from the ID line -> fix from the vertex sums later
    return tot


def main():
    n = int(sys.argv[1])
    i = sys.argv.index('--exact') if '--exact' in sys.argv else sys.argv.index('--kf')
    files = sys.argv[2:i]
    plug = load_plug(n, files)
    ident = plug.pop(((1,) * n, None), Fraction(0))
    print(f'n={n}: total weight {sum(plug.values()) + ident} (= A000944(n)), cells {len(plug)} + identity')
    if sys.argv[i] == '--exact':
        ref = defaultdict(Fraction)
        for nn, vp, eo, w in json.load(open(sys.argv[i + 1])):
            if nn == n:
                ref[(tuple(vp), tuple(tuple(e) for e in eo))] += Fraction(w)
        rid = sum(w for (vp, eo), w in ref.items() if all(x == 1 for x in vp))
        bad = sum(1 for k in set(plug) | {k for k in ref if not all(x == 1 for x in k[0])}
                  if plug.get(k, 0) != ref.get(k, 0))
        bad += ident != rid
        print('exact cells:', 'ALL MATCH' if not bad else f'{bad} MISMATCHES')
        return
    N, npr = int(sys.argv[i + 1]), int(sys.argv[i + 2])
    import tower_p as T
    allbad = 0
    for q in T.PRIMES24[:npr]:
        ref = defaultdict(int)
        rid = 0
        for nn, vp, eo, w in json.load(open(f'../kf_{N}_{q}.json')):
            if nn != n:
                continue
            if all(x == 1 for x in vp):
                rid = (rid + w) % q
            else:
                k = (tuple(vp), tuple(sorted(tuple(e) for e in eo)))
                ref[k] = (ref[k] + w) % q
        fm = lambda f: f.numerator % q * pow(f.denominator % q, q - 2, q) % q
        keys = set(plug) | {k for k, v in ref.items() if v}
        bad = [k for k in keys if fm(plug.get(k, Fraction(0))) != ref.get(k, 0) % q]
        bad += ['identity'] if fm(ident) != rid else []
        allbad += len(bad)
        print(q, f'{len(keys) - len(bad) + (fm(ident) == rid)}/{len(keys) + 1} cells match', bad[:3])
    print('ALL MATCH' if not allbad else f'{allbad} MISMATCHES')


if __name__ == '__main__':
    main()
