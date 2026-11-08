"""Collect the clearly symmetric panels (score >= 0.9, size >= 30), symmetrize each by majority vote over
the 8 symmetries of the square (ties -> open), and save them as panels/<page>_<k>.npy."""
import json, os, numpy as np
os.makedirs('panels', exist_ok=True)
res = json.load(open('panels.json')); n = 0
for f, ps in res.items():
    G = np.load(f)
    for k, p in enumerate(ps):
        if p['score'] < 0.9 or p['s'] < 30: continue
        A = G[p['r']:p['r'] + p['s'], p['c']:p['c'] + p['s']].astype(int)
        imgs = [np.rot90(A, q) for q in range(4)]; imgs += [x.T for x in imgs]
        S = (np.mean(imgs, 0) > 0.5).astype(np.uint8)
        np.save(f'panels/{f[:4]}_{k}.npy', S); n += 1
        print(f[:4], k, p['s'], p['score'], 'changed', int((S != A).sum()))
print(n, 'panels')
