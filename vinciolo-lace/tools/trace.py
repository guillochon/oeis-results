"""Trace a Vinciolo grid-pattern page photo into a 0/1 array (1 = printed/filled square).

Dark squares are found as connected blobs, then placed on an integer lattice by walking between neighboring
blobs (steps rounded to the local cell pitch), which follows page curvature. Cells with no blob are open.
    python trace.py page.jpg out.npy [threshold]
"""
import sys
from collections import deque
import numpy as np
from PIL import Image
from scipy import ndimage as nd
from scipy.spatial import cKDTree


def trace(path, thr=None):
    im = np.asarray(Image.open(path).convert('L')).astype(float)
    if thr is None:                       # midway between ink and paper levels
        lo, hi = np.percentile(im, [5, 90])
        thr = (lo + hi) / 2 - 0.15 * (hi - lo)
    mask = im < thr
    lab, n = nd.label(mask)
    sz = np.bincount(lab.ravel())[1:]
    cen = np.array(nd.center_of_mass(mask, lab, np.arange(1, n + 1)))
    med = np.median(sz[sz > 10])
    keep = (sz > 0.3 * med) & (sz < 2.2 * med)
    c = cen[keep]
    t = cKDTree(c)
    d, _ = t.query(c, k=2)
    pitch0 = np.median(d[:, 1])
    pairs = t.query_pairs(3.3 * pitch0, output_type='ndarray')
    v = c[pairs[:, 1]] - c[pairs[:, 0]]
    # refine per-axis pitch from near-axis unit steps
    dd = np.hypot(*v.T)
    unit = (dd > 0.7 * pitch0) & (dd < 1.3 * pitch0)
    py = np.median(np.abs(v[unit & (np.abs(v[:, 1]) < 0.3 * pitch0), 0]))
    px = np.median(np.abs(v[unit & (np.abs(v[:, 0]) < 0.3 * pitch0), 1]))
    s = np.array([py, px])
    q = v / s
    r = np.rint(q)
    good = (np.abs(q - r).max(1) < 0.2) & (np.abs(r).max(1) <= 3) & (r != 0).any(1)
    adj = [[] for _ in c]
    for (a, b), rr in zip(pairs[good], r[good].astype(int)):
        adj[a].append((b, rr))
        adj[b].append((a, -rr))
    coord = {}
    conflicts = 0
    comps = []
    for s0 in range(len(c)):
        if s0 in coord:
            continue
        coord[s0] = np.array([0, 0])
        comp = [s0]
        dq = deque([s0])
        while dq:
            a = dq.popleft()
            for b, rr in adj[a]:
                cb = coord[a] + rr
                if b in coord:
                    conflicts += (coord[b] != cb).any()
                else:
                    coord[b] = cb
                    comp.append(b)
                    dq.append(b)
        comps.append(comp)
    main = max(comps, key=len)
    cs = np.array([coord[i] for i in main])
    mn = cs.min(0)
    G = np.zeros(tuple(cs.max(0) - mn + 1), np.uint8)
    G[tuple((cs - mn).T)] = 1
    centers = c[main]
    info = dict(blobs=len(c), placed=len(main), conflicts=int(conflicts), pitch=s.tolist(), shape=G.shape)
    return G, info, (centers, cs - mn)


if __name__ == '__main__':
    G, info, _ = trace(sys.argv[1], float(sys.argv[3]) if len(sys.argv) > 3 else None)
    np.save(sys.argv[2], G)
    print(info)
