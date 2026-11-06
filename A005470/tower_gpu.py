"""GPU (CuPy) version of the tower ring and of the 3-connected substitution (tower_p.py equations).

PSG: the cycle-index series ring on the GPU, same interface as tower_p.PS (so tower_p's network,
block and univariate code runs unchanged on it).  Products use the precomputed pair list (I, J)
grouped by output index K: one thread per output partition sums its pairs (a_i b_j + a_j b_i), with
Barrett reduction, no atomics.  A weight cap w restricts to the K-order prefix of outputs of weight <= w.

CoreG: the moment tables compiled for the GPU.  For a table with rows vp (vertex partitions) and
entries (key, c), the substitution is  sum_vp p^vp * sum_key c * prod(key).  Kernel 1 forms
T[vp, i] = sum_key c * prod(key)[i] for the partitions i of weight <= N - |vp| (a prefix); kernel 2
adds T[vp, i] into the partition i + vp (precomputed index map), skipping weights above the cap.

Usage: python tower_gpu.py N nprimes   (needs kf_N_<p>.json; same output as run_tower.py)
"""
import json
import sys
import time

import cupy as cp
import numpy as np

import tower_p as T

_SRC = r'''
__device__ __forceinline__ unsigned long long redp(unsigned long long x, unsigned long long p,
                                                   unsigned long long m) {
    unsigned long long q = __umul64hi(x, m);
    unsigned long long r = x - q * p;
    while (r >= p) r -= p;
    return r;
}

extern "C" __global__ void pmul(const long long* a, const long long* b, const int* I, const int* J,
    const long long* seg, const int* kout, long long* out, int nseg,
    unsigned long long p, unsigned long long m)
{
    int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= nseg) return;
    unsigned long long acc = 0;
    int cnt = 0;
    for (long long t = seg[s]; t < seg[s + 1]; t++) {
        int i = I[t], j = J[t];
        unsigned long long v = (unsigned long long)a[i] * (unsigned long long)b[j];
        if (i != j) v += (unsigned long long)a[j] * (unsigned long long)b[i];
        acc += v;
        if (++cnt == 4096) { acc = redp(acc, p, m); cnt = 0; }
    }
    out[kout[s]] = (long long)redp(acc, p, m);
}

extern "C" __global__ void rowsum(const long long* Pm, int P, const long long* rowptr,
    const int* ekey, const long long* ecoef, const long long* off, const int* rowof,
    long long total, long long* Tout, unsigned long long p, unsigned long long m)
{
    long long pos = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (pos >= total) return;
    int r = rowof[pos];
    long long i = pos - off[r];
    unsigned long long acc = 0;
    int cnt = 0;
    for (long long e = rowptr[r]; e < rowptr[r + 1]; e++) {
        acc += (unsigned long long)ecoef[e] * (unsigned long long)Pm[(long long)ekey[e] * P + i];
        if (++cnt == 4096) { acc = redp(acc, p, m); cnt = 0; }
    }
    Tout[pos] = (long long)redp(acc, p, m);
}

extern "C" __global__ void scatter(const long long* Tv, const int* dst, const int* dwt, int cap,
    long long total, unsigned long long* out)
{
    long long pos = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (pos >= total) return;
    if (dwt[pos] > cap) return;
    atomicAdd(out + dst[pos], (unsigned long long)Tv[pos]);
}
'''
_mod = cp.RawModule(code=_SRC)
_pmul = _mod.get_function('pmul')
_rowsum = _mod.get_function('rowsum')
_scatter = _mod.get_function('scatter')


def _grid(n, thr=256):
    return ((n + thr - 1) // thr,), (thr,)


class PSG:
    def __init__(self, N, p):
        cpu = T.PS(N, p)
        self.cpu = cpu
        self.N, self.p, self.P = N, p, cpu.P
        self.parts, self.idx = cpu.parts, cpu.idx
        self.wt_host = cpu.wt
        self.wt = cp.asarray(cpu.wt)
        self.m = np.uint64((1 << 64) // p)
        self.I = cp.asarray(cpu.I.astype(np.int32))
        self.J = cp.asarray(cpu.J.astype(np.int32))
        seg = np.r_[cpu.Kstart, len(cpu.K)].astype(np.int64)
        self.seg = cp.asarray(seg)
        self.kout = cp.asarray(cpu.K[cpu.Kstart].astype(np.int32))
        self.cap_starts = cpu.cap_starts
        self.pl = {k: (cp.asarray(s), cp.asarray(d)) for k, (s, d) in cpu.pl.items()}
        self.d1src, self.d1dst = cp.asarray(cpu.d1src), cp.asarray(cpu.d1dst)
        self.d1mul = cp.asarray(cpu.d1mul.astype(np.int64))
        self.set_cap(N)

    def set_cap(self, w):
        self.w = w
        self.nseg = int(self.cap_starts[w])

    def zero(self):
        return cp.zeros(self.P, dtype=cp.int64)

    def one(self):
        a = self.zero()
        a[0] = 1
        return a

    def mono(self, parts):
        a = self.zero()
        key = tuple(sorted(parts, reverse=True))
        if sum(key) <= self.N:
            a[self.idx[key]] = 1
        return a

    def mul(self, a, b):
        out = self.zero()
        if self.nseg:
            a = cp.ascontiguousarray(a % self.p) if a.dtype != cp.int64 else a
            _pmul(*_grid(self.nseg), (a, b, self.I, self.J, self.seg, self.kout, out, np.int32(self.nseg),
                                      np.uint64(self.p), self.m))
        return out

    def pleth(self, a, k):
        src, dst = self.pl[k]
        out = self.zero()
        out[dst] = a[src]
        return out

    def dp1(self, a):
        out = self.zero()
        out[self.d1dst] = a[self.d1src] * self.d1mul % self.p
        return out

    def inv(self, c):
        return pow(int(c) % self.p, self.p - 2, self.p)

    def exp0(self, f):
        r, term = self.one(), self.one()
        for j in range(1, self.N + 1):
            term = self.mul(term, f) * self.inv(j) % self.p
            if not bool(term.any()):
                break
            r = (r + term) % self.p
        return r

    def powr(self, a, e):
        r = self.one()
        for _ in range(e):
            r = self.mul(r, a)
        return r


class CoreG:
    def __init__(self, ps, R, Hid, Htau):
        self.ps = ps
        N, p = ps.N, ps.p
        tabs = [T.moment_table(x, N, p) for x in (R, Hid, Htau)]
        keys = sorted({k for tab in tabs for inner in tab.values() for k, c in inner.items() if c},
                      key=lambda k: (len(k), k))
        self.keys = keys
        kidx = {k: i for i, k in enumerate(keys)}
        nP = np.searchsorted(ps.wt_host, np.arange(N + 1), side='right')
        self.tabs = []
        for tab in tabs:
            vps = sorted(vp for vp, inner in tab.items() if any(inner.values()))
            rowptr, ekey, ecoef, off, rowof, dst, dwt = [0], [], [], [0], [], [], []
            for r, vp in enumerate(vps):
                for k, c in tab[vp].items():
                    if c:
                        ekey.append(kidx[k])
                        ecoef.append(c)
                rowptr.append(len(ekey))
                L = int(nP[N - sum(vp)])
                off.append(off[-1] + L)
                rowof.append(np.full(L, r, dtype=np.int32))
                d = np.array([ps.idx[tuple(sorted(ps.parts[i] + vp, reverse=True))] for i in range(L)],
                             dtype=np.int32)
                dst.append(d)
                dwt.append(ps.wt_host[d].astype(np.int32))
            total = off[-1]
            self.tabs.append(dict(
                rowptr=cp.asarray(np.array(rowptr, dtype=np.int64)),
                ekey=cp.asarray(np.array(ekey, dtype=np.int32)),
                ecoef=cp.asarray(np.array(ecoef, dtype=np.int64)),
                off=cp.asarray(np.array(off, dtype=np.int64)),
                rowof=cp.asarray(np.concatenate(rowof) if rowof else np.zeros(0, np.int32)),
                dst=cp.asarray(np.concatenate(dst) if dst else np.zeros(0, np.int32)),
                dwt=cp.asarray(np.concatenate(dwt) if dwt else np.zeros(0, np.int32)),
                total=total))

    def evaluate(self, Did, Dtau):
        ps, p = self.ps, self.ps.p
        one = ps.one()
        d = {0: (Did - one) % p, 1: (Dtau - one) % p}
        powcache, prodcache = {}, {}

        def dpow(kind, j):
            key = (kind, j)
            if key not in powcache:
                powcache[key] = one if j == 0 else ps.mul(dpow(kind, j - 1), ps.pleth(d[kind[1]], kind[0]))
            return powcache[key]

        def prod(key):
            if key not in prodcache:
                prodcache[key] = one if not key else ps.mul(prod(key[:-1]), dpow(*key[-1]))
            return prodcache[key]

        Pm = cp.stack([prod(k) for k in self.keys]) if self.keys else cp.zeros((1, ps.P), cp.int64)
        outs = []
        for t in self.tabs:
            total = t['total']
            acc = cp.zeros(ps.P, dtype=cp.uint64)
            if total:
                Tv = cp.empty(total, dtype=cp.int64)
                _rowsum(*_grid(total), (Pm, np.int32(ps.P), t['rowptr'], t['ekey'], t['ecoef'], t['off'],
                                        t['rowof'], np.int64(total), Tv, np.uint64(p), ps.m))
                _scatter(*_grid(total), (Tv, t['dst'], t['dwt'], np.int32(ps.w), np.int64(total), acc))
            outs.append((acc % np.uint64(p)).astype(cp.int64))
        return outs


# ---------------------------------------------------------------- univariate stage (vectorized, CPU)

def evalp_fast(ps, Z, f):
    """Z with p_j -> f(x^j), all partitions at once: products built by number of parts."""
    N, p = ps.N, ps.p
    parts = ps.parts
    P = len(parts)
    G = np.zeros((N + 1, N + 1), dtype=np.int64)     # G[j] = f(x^j)
    for j in range(1, N + 1):
        G[j] = T.upleth(np.asarray(f, dtype=np.int64) % p, j, N)
    V = np.zeros((P, N + 1), dtype=np.int64)
    V[0, 0] = 1
    bylen = {}
    for i, lam in enumerate(parts):
        if lam:
            bylen.setdefault(len(lam), []).append(i)
    for L in sorted(bylen):
        ids = np.array(bylen[L])
        par = np.array([ps.idx[parts[i][:-1]] for i in ids])
        last = np.array([parts[i][-1] for i in ids])
        A, Bm = V[par], G[last]
        out = np.zeros_like(A)
        for s in range(N + 1):
            if Bm[:, s].any():
                out[:, s:] = (out[:, s:] + A[:, :N + 1 - s] * Bm[:, s:s + 1]) % p
        V[ids] = out
    Z = np.asarray(Z, dtype=np.int64) % p
    res = np.zeros(N + 1, dtype=np.int64)
    for c0 in range(0, P, 4096):
        res = (res + (Z[c0:c0 + 4096, None] * V[c0:c0 + 4096] % p).sum(axis=0)) % p
    return res


def run_mod_gpu(N, p, rows, verbose=False):
    t0 = time.time()
    ps = PSG(N, p)
    R, Hid, Htau = T.rows_from_zk(rows, p)
    core = CoreG(ps, R, Hid, Htau)
    if verbose:
        print(f'  setup {time.time() - t0:.0f}s', flush=True)
    B, K2 = T.blocks(ps, core, verbose)
    cp.cuda.Device().synchronize()
    if verbose:
        print(f'  blocks {time.time() - t0:.0f}s', flush=True)
    B, K2 = cp.asnumpy(B), cp.asnumpy(K2)
    T.evalp = evalp_fast                     # the univariate stage reads ps.parts/idx only
    cpu = ps.cpu
    one_x = np.zeros(N + 1, dtype=np.int64)
    one_x[1] = 1
    bc = evalp_fast(cpu, (B - K2) % p, one_x)
    C, Cu = T.connected(cpu, B)
    G = T.euler(Cu, N, p)
    if verbose:
        print(f'  done {time.time() - t0:.0f}s', flush=True)
    return bc, Cu, G


def main():
    N, npr = int(sys.argv[1]), int(sys.argv[2])
    primes = T.PRIMES24[:npr]
    res = []
    for q in primes:
        rows = json.load(open(f'kf_{N}_{q}.json'))
        print('prime', q, flush=True)
        res.append(run_mod_gpu(N, q, rows, verbose=True))
    import run_tower
    run_tower.report(N, primes, res)


if __name__ == '__main__':
    main()
