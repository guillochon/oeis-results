"""Search for non-local / growth rules that all 19 Vinciolo panels obey, and how often density-matched random
no-pinhole designs obey them.  Arrays: 1 = darned (printed), 0 = open; outside the panel counts as darned."""
import sys, numpy as np
from scipy import ndimage as nd
from growth_load import load
S8 = np.ones((3, 3), int); S4 = nd.generate_binary_structure(2, 1)
P = {k: np.pad(A, 1, constant_values=1) for k, A in load().items()}
def nulls(k, n):
    out = []
    for l in open(f'null/{k}.txt'):
        m = n - 2; B = np.array([int(c) for c in l.strip()]).reshape(m, m)
        out.append(np.pad(B, 2, constant_values=1))
    return out
N = {k: nulls(k, A.shape[0] - 2) for k, A in P.items()}

def radius(A, kind):
    n = A.shape[0]; c = (n - 1) / 2; i, j = np.indices(A.shape); di, dj = abs(i - c), abs(j - c)
    return {'cheb': np.maximum(di, dj), 'manh': di + dj, 'eucl': np.hypot(di, dj)}[kind]

def growth_viol(A, color, kind, outward, S):
    """cells of `color` lacking a same-color neighbor (in S) strictly closer to (outward=True: the center;
    False: the outside) -- except the extreme cells of their own component (the component's seed ring)."""
    R = radius(A, kind); M = A == color
    lab, k = nd.label(M, S); v = 0
    offs = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a or b) and S[a + 1, b + 1]]
    n = A.shape[0]
    for c in range(1, k + 1):
        C = lab == c; ext = R[C].min() if outward else R[C].max()
        for i, j in zip(*np.nonzero(C)):
            if np.isclose(R[i, j], ext): continue
            ok = False
            for a, b in offs:
                u, w = i + a, j + b
                if 0 <= u < n and 0 <= w < n and M[u, w] and (R[u, w] < R[i, j] - 1e-9 if outward else R[u, w] > R[i, j] + 1e-9):
                    ok = True; break
            v += not ok
    return v

def comp_sizes(A, color, S, touching_frame=True):
    lab, k = nd.label(A == color, S); sz = nd.sum(np.ones_like(A), lab, range(1, k + 1)).astype(int)
    if not touching_frame:
        edge = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])) - {0}
        sz = [s for c, s in enumerate(sz, 1) if c not in edge]
    return list(sz)

def eight_not_four(A, color):
    """components joined only through a diagonal pinch: #8-components < #4-components"""
    return nd.label(A == color, S4)[1] - nd.label(A == color, S8)[1]

def holes_in_open(A):
    """darned islands (not touching the outside) -- each sits inside an open region"""
    return len(comp_sizes(A, 1, S8, False))

def open_holes(A):
    """open regions whose 4-connected body encloses darned cells: Euler-type count via filling"""
    v = 0
    lab, k = nd.label(A == 0, S8)
    for c in range(1, k + 1):
        C = lab == c; F = nd.binary_fill_holes(C, S4)
        v += (F & ~C).any()
    return v

def main():
    tests = {}
    for col, cn in ((0, 'open'), (1, 'darned')):
        for kind in ('cheb', 'manh', 'eucl'):
            for out in (True, False):
                for S, sn in ((S8, '8'), (S4, '4')):
                    tests[f'{cn} grows {"outward from center" if out else "inward from edge"} ({kind}, {sn}-nbr)'] = \
                        (lambda A, col=col, kind=kind, out=out, S=S: growth_viol(A, col, kind, out, S), 0)
    tests['open 8-regions joined only diagonally'] = (lambda A: eight_not_four(A, 0), 0)
    tests['darned 8-regions joined only diagonally'] = (lambda A: eight_not_four(A, 1), 0)
    tests['open regions enclosing darned cells'] = (open_holes, 0)
    for k in (2, 3, 4, 5, 6, 8):
        tests[f'open 8-regions smaller than {k}'] = (lambda A, k=k: sum(s < k for s in comp_sizes(A, 0, S8)), 0)
        tests[f'open 4-regions smaller than {k}'] = (lambda A, k=k: sum(s < k for s in comp_sizes(A, 0, S4)), 0)
        tests[f'floating darned 8-islands of size 2..{k}'] = (lambda A, k=k: sum(2 <= s <= k for s in comp_sizes(A, 1, S8, False)), 0)

    rows = []
    for name, (f, target) in tests.items():
        pv = [f(A) for A in P.values()]
        held = sum(x == target for x in pv)
        nv = [f(B) == target for k in P for B in N[k]]
        rows.append((held, np.mean(nv), name, pv))
    rows.sort(key=lambda r: (-r[0], r[1]))
    for held, nr, name, pv in rows:
        print(f'{held:2d}/19  random {nr:6.1%}  {name:55s} {pv if held >= 15 else ""}')

if __name__ == '__main__':
    main()
