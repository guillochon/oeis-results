"""Independent check of the k = 4 computation (a(m+4) from L = J(M(m))) modulo one prime.

Second implementation, different from rust/src/main.rs in three ways: numpy on whole blocks of
columns, the dense pair table B(I_14, I_23) (all pairs, not only comparable ones), and the final
contraction sum_{e,f} B(e,f) B(f*,e*). Tables go to disk in blocks (about 16 GB for m = 7).

Run: python attacks/A132183/verify_dense.py M TMPDIR [BLOCK]
Prints a(M+4) mod 2^31-1 and the self-dual count mod 2^31-1.
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proto  # noqa: E402

P = 2147483647
m = int(sys.argv[1])
tmp = sys.argv[2]
BLK = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
os.makedirs(tmp, exist_ok=True)
t0 = time.time()

try:  # below-normal priority, and no EcoQoS throttling
    import psutil
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
except Exception:
    pass

L = proto.Lat(m)
N = L.N
print(f"m={m}: {N} ideals, {time.time()-t0:.1f} s", flush=True)
mask64 = (1 << 64) - 1
words = max(1, (1 << m) // 64)
idw = np.array([[(I >> (64 * w)) & mask64 for w in range(words)] for I in L.ideals], dtype=np.uint64)
sig, dual = L.sig, L.dual


def contains(rows, cols):
    """bool[len(rows), len(cols)]: ideal rows[i] contains ideal cols[j]."""
    out = np.ones((len(rows), len(cols)), dtype=bool)
    for w in range(words):
        out &= (idw[cols, w][None, :] & ~idw[rows, w][:, None]) == 0
    return out


def up(M):
    for p in reversed(range(L.n)):
        e = L.edges[p]
        if len(e):
            M[e[:, 1]] = (M[e[:, 1]] + M[e[:, 0]]) % P
    return M


def down(M):
    for p in range(L.n):
        e = L.edges[p]
        if len(e):
            M[e[:, 0]] = (M[e[:, 0]] + M[e[:, 1]]) % P
    return M


u = L.up(proto.ones(N))[sig]
v = L.up(u)
vmod = np.array([int(x) % P for x in v], dtype=np.uint64)
allidx = np.arange(N)
order = np.argsort(sig, kind="stable")
starts = np.flatnonzero(np.r_[True, sig[order][1:] != sig[order][:-1]])
sigvals = sig[order][starts]

# step A: T4[e, d] for blocks of d
dblocks = [np.arange(a, min(a + BLK, N)) for a in range(0, N, BLK)]
sd = 0
for bi, D in enumerate(dblocks):
    CD = contains(allidx, D)                                  # I >= d
    cnt = np.zeros((N, len(D)), dtype=np.uint32)
    cnt[sigvals] = np.add.reduceat(CD[order].astype(np.uint32), starts, axis=0)
    n = down(cnt)                                             # #{I_12 >= d : sigma(I_12) <= I_2}
    x = ((n.astype(np.uint64) * vmod[:, None]) % P).astype(np.uint32)
    T2 = up(x)
    T2 *= contains(allidx, sig[D])                            # I_3 >= sigma(d)
    T3p = up(up(T2))
    T4 = T3p[sig] * contains(D, allidx).T                     # e <= d
    sd = (sd + int((T4 * contains(allidx, dual[D]).astype(np.uint64)).sum() % P)) % P   # e* <= d
    np.save(os.path.join(tmp, f"T4_{bi}.npy"), T4)
    print(f"A block {bi+1}/{len(dblocks)} {time.time()-t0:.0f} s", flush=True)
print("selfdual mod p =", sd, flush=True)

# order in which the dual map preserves blocks
key = np.minimum(allidx, dual)
perm = np.lexsort((allidx, key))
posperm = np.empty(N, dtype=np.int64)
posperm[perm] = allidx
cuts = [0]
while cuts[-1] < N:
    c = min(cuts[-1] + BLK, N)
    while c < N and key[perm[c]] == key[perm[c - 1]]:
        c += 1
    cuts.append(c)
eblocks = [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]
T4f = [np.load(os.path.join(tmp, f"T4_{bi}.npy"), mmap_mode="r") for bi in range(len(dblocks))]

# step B: B[e, f] = sum over d >= f of T4[e, d]
for bi, (a, b) in enumerate(eblocks):
    E = perm[a:b]
    Es = np.sort(E)
    M = np.concatenate([np.asarray(f[Es, :]) for f in T4f], axis=1)   # [e, d]
    Bt = up(np.ascontiguousarray(M.T))                                # [f, e]
    back = np.searchsorted(Es, E)
    np.save(os.path.join(tmp, f"B_{bi}.npy"), np.ascontiguousarray(Bt.T[back][:, perm]))
    print(f"B block {bi+1}/{len(eblocks)} {time.time()-t0:.0f} s", flush=True)
del T4f
Bf = [np.load(os.path.join(tmp, f"B_{bi}.npy"), mmap_mode="r") for bi in range(len(eblocks))]
dloc = [posperm[dual[perm[a:b]]] - a for a, b in eblocks]
for (a, b), dl in zip(eblocks, dloc):
    assert dl.min() == 0 and dl.max() == b - a - 1

# final: sum over e, f of B[e, f] B[f*, e*]
acc = 0
for i, (a, b) in enumerate(eblocks):
    for j, (c, d) in enumerate(eblocks):
        X = np.asarray(Bf[i][:, c:d]).astype(np.uint64)
        Y = np.asarray(Bf[j][:, a:b])
        W = Y[dloc[j]][:, dloc[i]].T.astype(np.uint64)
        acc = (acc + int(((X * W) % P).sum())) % P
    print(f"final block {i+1}/{len(eblocks)} {time.time()-t0:.0f} s", flush=True)
print(f"a({m+4}) mod {P} = {acc}")
print(f"selfdual({m+4}) mod {P} = {sd}")
