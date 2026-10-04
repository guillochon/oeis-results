"""a(n) mod m for the stacking sequence (distinct sizes) is eventually periodic for every m: the slot
recurrence only multiplies by the slot counts (s, f, e, t), so it runs on the m^4 residue classes of the state.
This finds the exact preperiod and period by iterating the collapsed state until it repeats.

    python periodicity.py 2 3 4
"""
import sys
from collections import defaultdict


def step_mod(vec, m):
    def place(v, kind):
        out = defaultdict(int)
        for (s, f, e, t), w in v.items():
            if kind == 'D':
                moves = [(1, (s+1, f, e+1, t)), (s, (s, f, e+1, t)), (f, (s+1, f-1, e+1, t)), (e, (s, f, e, t))]
            elif kind == 'B':
                moves = [(1, (s, f+1, e, t)), (s, (s-1, f+1, e, t)), (f+e, (s, f, e, t))]
            else:
                moves = [(1, (s+1, f, e, t+1)), (s, (s, f, e, t+1)), (t, (s, f, e, t))]
            for mult, st in moves:
                k = tuple(x % m for x in st)
                out[k] = (out[k] + w * mult) % m
        return out
    d = place(vec, 'D')
    for k, w in place(place(vec, 'B'), 'T').items():
        d[k] = (d[k] + w) % m
    return {k: w for k, w in d.items() if w}


def period_mod(m, limit):
    v = {(0, 0, 0, 0): 1}
    seen, seq = {}, []
    for n in range(limit):
        key = tuple(sorted(v.items()))
        if key in seen:
            return seen[key], n - seen[key], seq
        seen[key] = n
        seq.append(sum(v.values()) % m)
        v = step_mod(v, m)
    return None, None, seq


if __name__ == '__main__':
    a = [int(l.split()[1]) for l in open('bfile_stacking.txt')]
    for m in map(int, sys.argv[1:]):
        pre, per, seq = period_mod(m, 20000)
        ok = all(seq[n] == a[n] % m for n in range(min(len(seq), len(a))))
        print(f'mod {m}: preperiod {pre}, period {per} (agrees with the b-file: {ok})', flush=True)
        if per and per <= 60:
            print('  one period:', ''.join(map(str, seq[pre:pre + per])), flush=True)
