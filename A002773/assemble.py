"""Assemble the matroid family at n = 9, 10 and check it.

Inputs (this directory):
  orbits_hyp_out.txt  from orbits_hyp.py on the Mayhew-Royle catalogue (n <= 9): simple counts,
                      Aut-orbit counts of parallel extensions, labeled counts by rank.
  jo_result.json      Jo 2026 (github.com/ainta/matroid-count v2.0.0): all matroids on 10 elements by rank,
                      unlabeled and labeled.
  jo_cmp.json         Jo's comparison with the catalogue: unlabeled counts by (n, rank), n <= 9.

Identities (U = unlabeled, L = labeled; r = rank):
  U loopless(n, r) = U all(n, r) - U all(n-1, r)              (delete one loop)
  U loopless(n, r) = sum_{k <= n} sum_{simple M, k points, rank r} #(Aut(M)-orbits on compositions of n into k parts)
  L all(n, r)      = sum_j C(n, j) L loopless(n-j, r)         (choose the loops)
  L loopless(n, r) = sum_k S2(n, k) L simple(k, r)            (partition into parallel classes)
Writes results.json and prints everything.
"""
import json
from math import comb

from sympy.functions.combinatorial.numbers import stirling

out = open('orbits_hyp_out.txt').read().splitlines()


def rows(prefix, key):
    d = {}
    for line in out:
        if line.startswith(prefix):
            n = int(line.split(key)[1].split(':')[0])
            d[n] = list(map(int, line.split(':')[1].split('sum')[0].split()))
    return d


simple_cat = rows('simple', 'k=')          # unlabeled simple, by rank 0..9, k <= 9
orbit = rows('loopless', 's=')             # unlabeled loopless built from simple k <= 9, by rank 0..10
lab_cat = rows('labeled all', 'n=')        # labeled all matroids, n <= 9, by rank

cmp_ = json.load(open('jo_cmp.json'))['counts_by_size_rank']
jo = json.load(open('jo_result.json'))
N = 10
UA = {n: [cmp_.get('%d,%d' % (n, r), 0) for r in range(n + 1)] for n in range(10)}
UA[10] = [int(jo['by_rank'][str(r)]['unlabeled']) for r in range(11)]
LA = {n: lab_cat[n] for n in range(10)}
LA[10] = [int(jo['by_rank'][str(r)]['labeled']) for r in range(11)]
assert sum(UA[10]) == int(jo['unlabeled']) and sum(LA[10]) == int(jo['labeled'])

# --- unlabeled ---
U_loopless = {n: [UA[n][r] - (UA[n - 1][r] if 0 < n and r < n else 0) for r in range(n + 1)] for n in range(N + 1)}
U_simple = {}
for n in range(N + 1):
    if n <= 9:
        assert U_loopless[n] == orbit[n][:n + 1], n            # catalogue vs Burnside-free orbit count
        U_simple[n] = [simple_cat[n][r] if r < 10 else 0 for r in range(n + 1)]
    else:
        U_simple[n] = [U_loopless[n][r] - orbit[n][r] for r in range(n + 1)]
        assert min(U_simple[n]) >= 0

# --- labeled ---
L_loopless, L_simple = {}, {}
for n in range(N + 1):
    L_loopless[n] = [LA[n][r] - sum(comb(n, j) * (L_loopless[n - j][r] if r <= n - j else 0) for j in range(1, n + 1))
                     for r in range(n + 1)]
for n in range(N + 1):
    L_simple[n] = [L_loopless[n][r] - sum(int(stirling(n, k)) * (L_simple[k][r] if r <= k else 0) for k in range(n))
                   for r in range(n + 1)]
    assert min(L_simple[n]) >= 0 and min(L_loopless[n]) >= 0
# labeled simple at n <= 9 must also equal sum over simple M of n!/|Aut| -- implied; cross-check with Jo's catalogue
# labeled total and with the OEIS rows below.
jo_lab = json.load(open('jo_cmp.json')).get('labeled_by_size_rank')
if jo_lab:
    for n in range(10):
        assert LA[n] == [jo_lab['%d,%d' % (n, r)] for r in range(n + 1)], n

# --- checks against OEIS DATA (snapshot 2026-10-23) ---
oeis = {
    'A055545': ([sum(UA[n]) for n in range(11)], '1,2,4,8,17,38,98,306,1724,383172,3232000741644'),
    'A058718': ([sum(U_loopless[n]) for n in range(10)], '1,1,2,4,9,21,60,208,1418,381448'),
    'A002773': ([sum(U_simple[n]) for n in range(10)], '1,1,1,2,4,9,26,101,950,376467'),
    'A058731': ([U_simple[n][3] if n >= 3 else 0 for n in range(11)], '0,0,0,1,2,4,9,23,68,383,5249'),
    'A058719': ([U_loopless[n][3] for n in range(3, 11)], '1,3,9,25,70,217,950,8762'),
    'A058693': ([UA[n][3] for n in range(3, 11)], '1,4,13,38,108,325,1275,10037'),
    'A336704': ([UA[n][4] for n in range(4, 11)], '1,5,23,108,940,190214,4886380924'),
    'A058673': ([sum(LA[n]) for n in range(9)], '1,2,5,16,68,406,3807,75164,10607540'),
    'A058712': ([sum(L_loopless[n]) for n in range(9)], '1,1,2,6,27,185,2135,55129,10094077'),
    'A058721': ([sum(L_simple[n]) for n in range(2, 9)], '1,2,7,49,733,29760,9000402'),
    'A058722': ([L_simple[n][4] for n in range(4, 9)], '1,16,337,18700,7642631'),
    'A058715': ([L_loopless[n][3] for n in range(3, 11)], '1,11,106,1232,22172,803583,70820187,16122092568'),
    'A058681': ([LA[n][2] for n in range(2, 11)], '1,7,36,171,813,4012,20891,115463,677546'),
    'A100728': ([L_simple[n][n - 2] for n in range(4, 11)], '1,31,337,2570,16865,104858,650761'),
    'A058687': ([LA[n][3] for n in range(3, 9)], '1,15,171,2053,33442,1022217'),
    'A058669': ([x for n in range(9) for x in LA[n]] + LA[9][:2],
                '1,1,1,1,3,1,1,7,7,1,1,15,36,15,1,1,31,171,171,31,1,1,63,813,2053,813,63,1,1,127,4012,33442,33442,'
                '4012,127,1,1,255,20891,1022217,8520812,1022217,20891,255,1,1,511'),
    'A058710': ([x for n in range(9) for x in L_loopless[n]],
                '1,0,1,0,1,1,0,1,4,1,0,1,14,11,1,0,1,51,106,26,1,0,1,202,1232,642,57,1,0,1,876,22172,28367,3592,120,1,'
                '0,1,4139,803583,8274374,991829,19903,247,1'),
    'A058720': ([x for n in range(2, 9) for x in L_simple[n][2:]],
                '1,1,1,1,5,1,1,31,16,1,1,352,337,42,1,1,8389,18700,2570,99,1,1,433038,7642631,907647,16865,219,1'),
    'A058716': ([x for n in range(9) for x in U_loopless[n]],
                '1,0,1,0,1,1,0,1,2,1,0,1,4,3,1,0,1,6,9,4,1,0,1,10,25,18,5,1,0,1,14,70,85,31,6,1,0,1,21,217,832,288,51,7,1'),
    'A058730': ([x for n in range(2, 10) for x in U_simple[n][2:]] + U_simple[10][2:4],
                '1,1,1,1,2,1,1,4,3,1,1,9,11,4,1,1,23,49,22,5,1,1,68,617,217,40,6,1,1,383,185981,188936,1092,66,7,1,1,5249'),
    'A053534': ([x for n in range(10) for x in UA[n]],
                '1,1,1,1,2,1,1,3,3,1,1,4,7,4,1,1,5,13,13,5,1,1,6,23,38,23,6,1,1,7,37,108,108,37,7,1,1,8,58,325,940,325,'
                '58,8,1,1,9,87,1275,190214,190214,1275,87,9,1'),
}
for a, (mine, data) in oeis.items():
    assert mine == list(map(int, data.split(','))), (a, mine)
    print('%s: all %d known terms reproduced' % (a, len(mine)))
assert U_simple[10][4] == 4884573865  # A058730 row 10 (Matsumoto et al. 2012), = A058733(10)

res = {
    'unlabeled_all': UA, 'unlabeled_loopless': U_loopless, 'unlabeled_simple': U_simple,
    'labeled_all': LA, 'labeled_loopless': L_loopless, 'labeled_simple': L_simple,
}
json.dump({k: {str(n): v for n, v in d.items()} for k, d in res.items()}, open('results.json', 'w'), indent=1)
for k, d in res.items():
    for n in (9, 10):
        print('%-19s n=%2d total %d  by rank %s' % (k, n, sum(d[n]), ','.join(map(str, d[n]))))
print('A002773(10) =', sum(U_simple[10]))
