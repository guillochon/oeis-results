"""Independent check of a sparse witness: chi(q)=1 except E; chi(n) by trial division per chunk.
usage: python ls/verify_chunked.py k X p:v ...   -> reports first r<X with chi(r)=chi(r+1)=0"""
import sys, numpy as np
k, X = int(sys.argv[1]), int(sys.argv[2])
E = {int(a): int(b) for a, b in (s.split(":") for s in sys.argv[3:])}
chi = lambda q: E.get(q, 1) % k
R = int(X ** 0.5) + 2
s = np.ones(R + 1, bool); s[:2] = False
for i in range(2, int(R ** .5) + 1):
    if s[i]: s[i*i::i] = False
P = np.nonzero(s)[0]
C = 10_000_000; prev = False; first = None
for a in range(1, X + 1, C):
    b = min(a + C, X + 1)
    n = np.arange(a, b, dtype=np.int64); rem = n.copy(); c = np.zeros(b - a, np.int64)
    for p in P:
        p = int(p)
        if p * p > b and p > R: break
        st = (-a) % p
        idx = np.arange(st, b - a, p)
        while idx.size:
            rem[idx] //= p; c[idx] += chi(p)
            idx = idx[rem[idx] % p == 0]
    big = rem > 1
    for q in np.unique(rem[big]):
        if int(q) in E: c[rem == q] += chi(int(q)) - 1
    c[big] += 1
    z = (c % k) == 0
    if prev and z[0]: first = a - 1; break
    w = np.nonzero(z[:-1] & z[1:])[0]
    if w.size: first = a + int(w[0]); break
    prev = bool(z[-1]); print(f"ok to {b-1}", flush=True)
print("first zero pair r =", first)
