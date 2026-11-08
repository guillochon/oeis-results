"""Find square windows of a traced page that are (nearly) symmetric under all 8 symmetries of the square.
For every size s and position, score = mean agreement with its 7 nontrivial D4 images. Keep the best
non-overlapping windows with score >= MIN and s >= SMIN, preferring large ones."""
import sys, glob, json
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv
MIN, SMIN = 0.83, 21
def scan(G):
    G = G.astype(np.int8); out = []
    for s in range(min(G.shape), SMIN - 1, -1):
        W = swv(G, (s, s))                        # (H-s+1, W-s+1, s, s)
        imgs = [np.rot90(W, k, axes=(2, 3)) for k in range(4)]
        imgs += [np.swapaxes(x, 2, 3) for x in imgs]
        sc = np.mean([(W == B).mean(axis=(2, 3)) for B in imgs[1:]], axis=0)
        # require some structure: density between 0.15 and 0.85
        dens = W.mean(axis=(2, 3))
        sc[(dens < 0.15) | (dens > 0.85)] = 0
        # the window must sit on Vinciolo's solid frame: its outer ring at least 90% filled
        ring = (W[:, :, 0, :].sum(-1) + W[:, :, -1, :].sum(-1) + W[:, :, 1:-1, 0].sum(-1) + W[:, :, 1:-1, -1].sum(-1)) / (4 * s - 4)
        sc[ring < 0.9] = 0
        for i, j in zip(*np.nonzero(sc >= MIN)):
            out.append((s, float(sc[i, j]), int(i), int(j)))
    out.sort(key=lambda t: -t[0] * (t[1] - 0.8))   # big and symmetric
    taken = np.zeros(G.shape, bool); keep = []
    for s, sc, i, j in out:
        box = taken[i:i + s, j:j + s]
        if box.mean() > 0.2:                     # mostly inside a panel we already have
            continue
        keep.append(dict(s=s, score=round(sc, 3), r=i, c=j))
        taken[i:i + s, j:j + s] = True
    return keep
if __name__ == '__main__':
    res = {}
    for f in sys.argv[1:]:
        res[f] = scan(np.load(f))
        print(f, [(p['s'], p['score'], p['r'], p['c']) for p in res[f]])
    json.dump(res, open('panels.json', 'w'), indent=1)
