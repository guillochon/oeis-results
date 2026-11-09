"""Assemble T(n,k) by antidiagonals: T(1,1), T(1,2), T(2,1), T(1,3), T(2,2), T(3,1), ...
T(n,k) = number of n x n 0/1 arrays invariant under the 8 symmetries of the square in which every 8-connected
region of 1s has at least k cells.

Sources: k = 1 is 2^A008805(n-1) (every symmetric array); k = 2 is b_pinhole_d4.txt; k >= 3 is tools/table_runs.txt.
Cross-check: the brute-force table tools/minregion_counts.txt (n <= 12, k <= 8) must agree.
Writes b_minregion_table.txt (index m, T) for every complete antidiagonal, and prints the rows.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
T = {}
for l in open(os.path.join(ROOT, 'b_pinhole_d4.txt')):
    if l.strip() and not l.startswith('#'):
        n, v = l.split(); T[(int(n), 2)] = int(v)
def octant(n): h = (n + 1) // 2; return h * (h + 1) // 2
for n in range(1, 200): T[(n, 1)] = 2 ** octant(n)
for l in open(os.path.join(HERE, 'table_runs.txt')):
    n, k, v = map(int, l.split())
    if (n, k) in T and T[(n, k)] != v: raise SystemExit(f'conflict at {n},{k}')
    T[(n, k)] = v
for l in open(os.path.join(HERE, 'minregion_counts.txt')):
    if 'mode=0' in l:
        f = l.split(); n = int(f[0][2:])
        for k, v in enumerate(map(int, f[2:]), 1):
            assert T.get((n, k), v) == v, (n, k)
# monotonicity sanity check: T(n,k) is nonincreasing in k and T(n,k) = 1 for k > n^2
for (n, k), v in T.items():
    if (n, k + 1) in T: assert T[(n, k + 1)] <= v, (n, k)
    if k > n * n: assert v == 1, (n, k)
d = 2
seq = []
while all((d - k, k) in T for k in range(1, d)):
    seq += [T[(n, d - n)] for n in range(1, d)]; d += 1
last = d - 1
with open(os.path.join(ROOT, 'b_minregion_table.txt'), 'w', newline='\n') as f:
    for i, v in enumerate(seq, 1): f.write(f'{i} {v}\n')
print(f'complete antidiagonals n + k = 2..{last}: {len(seq)} terms')
print('next antidiagonal missing:', [(d - k, k) for k in range(1, d) if (d - k, k) not in T])
for n in range(1, 13):
    print(n, [T.get((n, k)) for k in range(1, 9)])
