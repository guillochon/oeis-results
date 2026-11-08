import glob, os, numpy as np
from scipy import ndimage as nd
P = {os.path.basename(f)[:-4]: np.load(f).astype(int) for f in sorted(glob.glob('panels/*.npy'))}
S4 = None; S8 = np.ones((3, 3))
def framed(A): return np.pad(A, 1, constant_values=1)
def comps(M, s): return nd.label(M, s)[1]
def lone(A, v, nb8=True):
    F = framed(A); n = A.shape[0]
    offs = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0) and (nb8 or a * b == 0)]
    same = sum((F[1 + a:n + 1 + a, 1 + b:n + 1 + b] == A) for a, b in offs)
    return int(((A == v) & (same == 0)).sum())
tests = {
 'open: one region (corner contact counts)': lambda A: comps(A == 0, S8),
 'open: one region (edges only)':           lambda A: comps(A == 0, S4),
 'filled+frame: one region (corners)':      lambda A: comps(framed(A) == 1, S8),
 'filled+frame: one region (edges only)':   lambda A: comps(framed(A) == 1, S4),
 'lone filled cells (8 nbrs)':              lambda A: lone(A, 1),
 'lone open cells (8 nbrs)':                lambda A: lone(A, 0),
 'lone filled cells (4 nbrs)':              lambda A: lone(A, 1, False),
 'lone open cells (4 nbrs)':                lambda A: lone(A, 0, False),
}
for name, f in tests.items():
    v = [f(A) for A in P.values()]
    target = 1 if 'one region' in name else 0
    print(f'{name:42s} holds in {sum(x == target for x in v):2d}/{len(v)}   values {sorted(v)}')
