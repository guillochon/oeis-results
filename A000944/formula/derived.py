"""b-files of the sequences that follow from the polyhedra tables U(v,f) and S(v,f) (v <= 40, f <= 76).

U = polyhedra (up to all isomorphisms), S = sensed (orientation-preserving isomorphisms only), achiral
= 2U - S. Duality U(v,f) = U(f,v), S(v,f) = S(f,v) reaches the cells with f <= 40, v <= 76. Every
b-file is checked term by term against the live OEIS terms (DATA and b-file, downloaded to LIVE) before
it is written.

    python derived.py LIVE_DIR      ->  ../derived/bNNNNNN.txt for each sequence, and a summary
"""
import json
import os
import sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'derived')
VX, FY = 40, 76
U = {tuple(map(int, k.split(','))): v for k, v in json.load(open(os.path.join(HERE, '../results/U_40_76.json'))).items()}
S = {tuple(map(int, k.split(','))): v for k, v in json.load(open(os.path.join(HERE, '../results/sensed_116.json'))).items()}


def polyhedral(v, f):
    """Cells (v vertices, f faces) that contain polyhedra: v, f >= 4, f <= 2v - 4, v <= 2f - 4."""
    return v >= 4 and f >= 4 and f <= 2 * v - 4 and v <= 2 * f - 4


def cell(T, v, f):
    if not polyhedral(v, f):
        return 0
    if v <= VX and f <= FY:
        return T.get((v, f), 0)
    if f <= VX and v <= FY:
        return T.get((f, v), 0)
    raise KeyError((v, f))


def u(v, f):
    return cell(U, v, f)


def s(v, f):
    return cell(S, v, f)


def achiral(v, f):
    return 2 * u(v, f) - s(v, f)


def by_edges(T, e):
    return sum(T(v, e + 2 - v) for v in range(1, e + 2))


def faces_row(T, n, ks):
    return [T(k, n) for k in ks]


# A-number -> (offset, list of terms, name of the range for the b-file line)
SEQS = {}
SEQS['A000109'] = (3, [1] + [u(n, 2 * n - 4) for n in range(4, 41)], 'n = 3..40')
SEQS['A002840'] = (6, [by_edges(u, e) for e in range(6, 80)], 'n = 6..79')
SEQS['A005645'] = (6, [by_edges(s, e) for e in range(6, 80)], 'n = 6..79')
SEQS['A006445'] = (1, [by_edges(achiral, e) for e in range(1, 80)], 'n = 1..79')
SEQS['A119501'] = (1, [sum(s(n, f) for f in range(4, 2 * n - 3)) for n in range(1, 41)], 'n = 1..40')
SEQS['A002856'] = (4, [u(n, n) for n in range(4, 41)], 'n = 4..40')
SEQS['A058789'] = (4, [u(n + 1, n) for n in range(4, 41)], 'n = 4..40')
SEQS['A342057'] = (4, [s(n, n) for n in range(4, 41)], 'n = 4..40')
# triangles, flattened
t = []
for n in range(4, 41):                                   # n faces, k = ceil(n/2) + 2 .. 2n - 4 vertices
    t += faces_row(u, n, range(-(-n // 2) + 2, 2 * n - 3))
SEQS['A058787'] = (4, t, None)
t = []
for n in range(4, 41):                                   # n faces, k = 4..2n-4 vertices
    t += faces_row(u, n, range(4, 2 * n - 3))
SEQS['A212438'] = (4, t, 'rows 4..40')
t = []
for n in range(4, 41):
    t += faces_row(s, n, range(4, 2 * n - 3))
SEQS['A239893'] = (4, t, 'rows 4..40')
t = []
for n in range(6, 80):                                   # n edges, k vertices with a polyhedron
    t += [u(k, n + 2 - k) for k in range(1, n + 2) if polyhedral(k, n + 2 - k)]
SEQS['A058788'] = (6, t, 'rows 6..79')
t = []
for n in range(1, 40):                                   # n nodes, k = 0..C(n,2) edges
    t += [u(n, k - n + 2) if k - n + 2 >= 1 else 0 for k in range(comb(n, 2) + 1)]
SEQS['A049337'] = (1, t, 'rows 1..39')


# Live terms we found to be wrong, with the reason (checked independently, see the README).
CORRECTIONS = {
    ('A005645', 34): 'digit typo 3323258053528 -> 3323268053528 (DATA and b-file); our U at 34 edges equals '
                     'A002840(34) from generation, and U = S/2 + W_ref + W_rr, so S = 3323258053528 would make it '
                     'differ by 5*10^6',
}

# A049337: the live DATA and b-file are the triangle of the name with its first three entries (zeros)
# missing, so every live term sits 3 places early (the example's row 4 puts K4 at k = 3 instead of 6).
SHIFTED = {'A049337': 3}


def live_terms(live, anum):
    """{n: a(n)} from the live b-file (if any) and DATA."""
    n = anum[1:]
    d = json.load(open(os.path.join(live, f'{anum}.json')))
    r = (d if isinstance(d, list) else d['results'])[0]
    off = int(r['offset'].split(',')[0])
    known = {off + i: int(x) for i, x in enumerate(r['data'].split(','))}
    bf = os.path.join(live, f'b{n}.txt')
    if os.path.exists(bf):
        for line in open(bf):
            p = line.split()
            if len(p) >= 2 and not line.startswith('#'):
                known[int(p[0])] = int(p[1])
    return known, r


def main():
    live = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    summary = {}
    for anum, (off, terms, rng) in SEQS.items():
        known, r = live_terms(live, anum)
        if anum in SHIFTED:
            sh = SHIFTED[anum]
            assert all(known[n] == terms[n - off + sh] for n in known), (anum, 'not a uniform shift')
            known = {n + sh: x for n, x in known.items()}
            print(f'{anum}: all {len(known)} live terms equal ours shifted by {sh} places (correction)')
        ours = {off + i: x for i, x in enumerate(terms)}
        overlap = [n for n in known if n in ours]
        bad = [n for n in overlap if known[n] != ours[n] and (anum, n) not in CORRECTIONS]
        fixed = [n for n in overlap if known[n] != ours[n] and (anum, n) in CORRECTIONS]
        assert not bad, (anum, bad[:5], [(known[n], ours[n]) for n in bad[:3]])
        assert set(known) <= set(ours), (anum, 'live terms beyond ours', sorted(set(known) - set(ours))[:5])
        n = anum[1:]
        with open(os.path.join(OUT, f'b{n}.txt'), 'w', newline='\n') as fh:
            for k in sorted(ours):
                fh.write(f'{k} {ours[k]}\n')
        summary[anum] = dict(offset=off, known=len(known), last_known=max(known), ours=len(ours),
                             last=max(ours), checked=len(overlap), corrected={n: [known[n], ours[n]] for n in fixed})
        print(f'{anum}: {len(overlap)} live terms checked, ' + (f'all equal except corrected {fixed}' if fixed else 'all equal') + f'; b-file n = {off}..{max(ours)} '
              f'({len(ours)} terms, {len(ours) - len(known)} new)')
    json.dump(summary, open(os.path.join(OUT, 'summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
