"""GPU port of the half-map grammar ('refl' mode of A000944/formula/disk.py) with three markers
(alpha on V labels, beta on X labels, gamma on V->V boundary edges), batched over marker points.

Same equations as disk.Grammar (arcs / virtuals / rarc / walk / root / cycle), rewritten for speed:
  * values are None (zero), Mono (x^i y^j times an optional per-batch scalar) or Dense
    ((B or 1) x (X+1) x (Y+1) int64 on the GPU) with a host-side lower bound `lo` on the x-degree,
    so zero tests never synchronize and products past the x truncation vanish structurally;
  * products Dense*Dense are an exact direct convolution kernel, Mono*Dense a shift-and-scale kernel,
    all reductions mod p by Barrett (no int64 division);
  * fixed-point pass t of the arc system is computed only for x-degree < t + 1: an arc of x-degree d
    is built from arcs of degree <= d - 1 (every nesting adds an interior vertex of weight >= x),
    so this loses nothing, and the last passes run at full size until the system is stable.
Validated against the CPU grammar (kformula.py check / identical rows).
"""
import cupy as cp
import numpy as np

V_, X_, P_ = 0, 1, 2
ANCHOR, ISV, LEXP = [1, 1, 0], [1, 0, 0], [1, 2, 2]

_SRC = r'''
__device__ __forceinline__ unsigned long long redp(unsigned long long x, unsigned long long p,
                                                   unsigned long long m) {
    unsigned long long q = __umul64hi(x, m);
    unsigned long long r = x - q * p;
    while (r >= p) r -= p;
    return r;
}

extern "C" __global__ void conv2(const long long* A, const long long* B, long long* C,
    int X, int Y, long long sa, long long sb, int loA, int loB, int cap,
    unsigned long long p, unsigned long long m)
{
    extern __shared__ unsigned long long sh[];
    int rows = cap;                       // rows 0..cap-1 are loaded (rows below lo are zero)
    unsigned long long* a = sh;
    unsigned long long* b = sh + rows * Y;
    const long long bt = blockIdx.x;
    const long long* Ag = A + bt * sa;
    const long long* Bg = B + bt * sb;
    for (int t = threadIdx.x; t < rows * Y; t += blockDim.x) {
        a[t] = (unsigned long long)Ag[t];
        b[t] = (unsigned long long)Bg[t];
    }
    __syncthreads();
    long long* Cg = C + bt * (long long)(X * Y);
    for (int o = threadIdx.x; o < X * Y; o += blockDim.x) {
        int i = o / Y, j = o % Y;
        unsigned long long acc = 0;
        if (i < cap && i >= loA + loB) {
            for (int k = loA; k <= i - loB; k++) {
                const unsigned long long* ar = a + k * Y;
                const unsigned long long* br = b + (i - k) * Y + j;
                for (int l = 0; l <= j; l++) acc += ar[l] * br[-l];
            }
            acc = redp(acc, p, m);
        }
        Cg[o] = (long long)acc;
    }
}

extern "C" __global__ void shiftmul(const long long* D, const long long* S, long long* C,
    int X, int Y, long long sd, int haveS, int mi, int mj, int cap, long long total,
    unsigned long long p, unsigned long long m)
{
    long long t = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (t >= total) return;
    long long bt = t / (X * Y);
    int o = (int)(t % (X * Y));
    int i = o / Y, j = o % Y;
    unsigned long long v = 0;
    if (i < cap && i >= mi && j >= mj) {
        v = (unsigned long long)D[bt * sd + (long long)(i - mi) * Y + (j - mj)];
        if (haveS) v = redp(v * (unsigned long long)S[bt], p, m);
    }
    C[t] = (long long)v;
}
'''
_mod = cp.RawModule(code=_SRC)
_conv2 = _mod.get_function('conv2')
_shift = _mod.get_function('shiftmul')

_add = cp.ElementwiseKernel('int64 a, int64 b, int64 p', 'int64 c',
                            'long long s = a + b; c = s >= p ? s - p : s;', 'gg_add')
_sub = cp.ElementwiseKernel('int64 a, int64 b, int64 p', 'int64 c',
                            'long long s = a - b; c = s < 0 ? s + p : s;', 'gg_sub')
_vmul = cp.ElementwiseKernel('int64 a, int64 b, int64 p', 'int64 c',
                             'c = (long long)(((unsigned long long)a * (unsigned long long)b) % (unsigned long long)p);',
                             'gg_vmul')


class Mono:
    __slots__ = ('i', 'j', 's')

    def __init__(self, i, j, s=None):
        self.i, self.j, self.s = i, j, s     # s: None (= 1) or cupy int64 vector (B,)


class Dense:
    __slots__ = ('a', 'lo')

    def __init__(self, a, lo):
        self.a, self.lo = a, lo              # a: (Bt, X+1, Y+1) int64, Bt = 1 or B


class Ring:
    def __init__(self, VX, FY, p):
        self.VX, self.FY, self.p = VX, FY, p
        self.X, self.Y = VX + 1, FY + 1
        self.cap = self.X
        self.m = np.uint64((1 << 64) // p)
        self.B = 1

    # ---------------------------------------------------------------- conversions
    def dense(self, v):
        if v is None or isinstance(v, Dense):
            return v
        if v.i >= self.X or v.j >= self.Y:
            return None
        Bt = 1 if v.s is None else v.s.shape[0]
        a = cp.zeros((Bt, self.X, self.Y), dtype=cp.int64)
        a[:, v.i, v.j] = 1 if v.s is None else v.s
        return Dense(a, v.i)

    def from_host(self, arr):
        arr = np.asarray(arr, dtype=np.int64) % self.p
        if not arr.any():
            return None
        lo = int(np.nonzero(arr.any(axis=tuple(range(arr.ndim - 1))))[0].min()) if arr.ndim == 2 else \
            int(np.nonzero(arr.reshape(-1, arr.shape[-2], arr.shape[-1]).any(axis=(0, 2)))[0].min())
        a = cp.asarray(arr.reshape((-1,) + arr.shape[-2:]))
        return Dense(a, lo)

    def to_host(self, v, B):
        if v is None:
            return np.zeros((B, self.X, self.Y), dtype=np.int64)
        d = self.dense(v)
        return np.broadcast_to(cp.asnumpy(d.a), (B, self.X, self.Y)).copy()

    # ---------------------------------------------------------------- arithmetic
    def add(self, u, v):
        if u is None:
            return v
        if v is None:
            return u
        if isinstance(u, Mono) and isinstance(v, Mono) and u.i == v.i and u.j == v.j:
            if u.s is None and v.s is None:
                return Mono(u.i, u.j, cp.full((1,), 2 % self.p, dtype=cp.int64))
            s1 = u.s if u.s is not None else cp.ones((1,), dtype=cp.int64)
            s2 = v.s if v.s is not None else cp.ones((1,), dtype=cp.int64)
            return Mono(u.i, u.j, _add(s1, s2, self.p))
        du, dv = self.dense(u), self.dense(v)
        if du is None:
            return dv
        if dv is None:
            return du
        return Dense(_add(du.a, dv.a, self.p), min(du.lo, dv.lo))

    def sub(self, u, v):
        if v is None:
            return u
        dv = self.dense(v)
        if dv is None:
            return u
        du = self.dense(u)
        if du is None:
            return Dense(_sub(cp.zeros_like(dv.a), dv.a, self.p), dv.lo)
        return Dense(_sub(du.a, dv.a, self.p), min(du.lo, dv.lo))

    def mul(self, u, v):
        if u is None or v is None:
            return None
        if isinstance(u, Mono) and isinstance(v, Mono):
            i, j = u.i + v.i, u.j + v.j
            if i >= min(self.cap, self.X) or j >= self.Y:
                return None
            if u.s is None:
                s = v.s
            elif v.s is None:
                s = u.s
            else:
                s = _vmul(u.s, v.s, self.p)
            return Mono(i, j, s)
        if isinstance(u, Mono):
            u, v = v, u
        if isinstance(v, Mono):
            return self._shift(u, v)
        return self._conv(u, v)

    def _shift(self, d, mo):
        lo = d.lo + mo.i
        cap = min(self.cap, self.X)
        if lo >= cap or mo.j >= self.Y:
            return None
        Bd = d.a.shape[0]
        Bs = 1 if mo.s is None else mo.s.shape[0]
        Bt = max(Bd, Bs)
        s = mo.s if mo.s is not None else cp.zeros((1,), dtype=cp.int64)
        if Bs == 1 and Bt > 1 and mo.s is not None:
            s = cp.broadcast_to(s, (Bt,)).copy()
        out = cp.empty((Bt, self.X, self.Y), dtype=cp.int64)
        total = Bt * self.X * self.Y
        thr = 256
        _shift(((total + thr - 1) // thr,), (thr,),
               (d.a, s, out, np.int32(self.X), np.int32(self.Y), np.int64(0 if Bd == 1 else self.X * self.Y),
                np.int32(0 if mo.s is None else 1), np.int32(mo.i), np.int32(mo.j), np.int32(cap),
                np.int64(total), np.uint64(self.p), self.m))
        return Dense(out, lo)

    def _conv(self, u, v):
        lo = u.lo + v.lo
        cap = min(self.cap, self.X)
        if lo >= cap:
            return None
        Ba, Bb = u.a.shape[0], v.a.shape[0]
        Bt = max(Ba, Bb)
        out = cp.empty((Bt, self.X, self.Y), dtype=cp.int64)
        shmem = 16 * cap * self.Y
        threads = min(1024, max(64, ((self.X * self.Y + 31) // 32) * 32))
        _conv2((Bt,), (threads,),
               (u.a, v.a, out, np.int32(self.X), np.int32(self.Y),
                np.int64(0 if Ba == 1 else self.X * self.Y), np.int64(0 if Bb == 1 else self.X * self.Y),
                np.int32(u.lo), np.int32(v.lo), np.int32(cap), np.uint64(self.p), self.m),
               shared_mem=shmem)
        return Dense(out, lo)

    def equal(self, u, v):
        du, dv = self.dense(u), self.dense(v)
        if du is None or dv is None:
            return du is None and dv is None or bool(not (du or dv).a.any())
        return bool(cp.array_equal(cp.broadcast_to(du.a, dv.a.shape) if du.a.shape[0] == 1 else du.a,
                                   cp.broadcast_to(dv.a, du.a.shape) if dv.a.shape[0] == 1 else dv.a))


def matmul(R, A, B):
    n, m, q = len(A), len(B), len(B[0])
    out = [[None] * q for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if A[i][k] is None:
                continue
            for j in range(q):
                if B[k][j] is not None:
                    out[i][j] = R.add(out[i][j], R.mul(A[i][k], B[k][j]))
    return out


def matadd(R, A, B):
    return [[R.add(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


class GGrammar:
    """cX: {d: host array (VX+1, FY+1)} of c-net weights by root face degree (disk.Grammar.cX)."""

    def __init__(self, VX, FY, p, cX):
        self.VX, self.FY, self.p = VX, FY, p
        self.R = R = Ring(VX, FY, p)
        self.cX = {d: R.from_host(a) for d, a in cX.items()}
        self.cX = {d: v for d, v in self.cX.items() if v is not None}
        self.dmax = max(cX) if cX else 3
        self.Fin = Mono(0, 2)
        self.Rarc0 = None
        self.L = 3

    def zeros(self):
        return [[None] * 3 for _ in range(3)]

    def omega(self, a, b, gam):
        j = 1 if ANCHOR[a] and not (ISV[a] and ISV[b]) else 0
        if ISV[a] and ISV[b] and gam is not None:
            return Mono(0, j, gam)
        return Mono(0, j)

    def lam(self, m):
        if m is None:
            return [None, None, Mono(2, 0)]
        a, b, _ = m
        return [Mono(1, 0, a), Mono(2, 0, b), Mono(2, 0)]

    # ------------------------------------------------------------------ arcs
    def arcs(self, m):
        """m: None (anchor-free; must come first, stores Rarc0) or (alpha, beta, gamma) cupy vectors."""
        R = self.R
        zero_m = m is None
        lam = self.lam(m)
        gam = None if zero_m else m[2]
        Real = [[R.mul(self.omega(a, b, gam), self.Fin) for b in range(3)] for a in range(3)]
        RealNoFace = [[self.omega(a, b, gam) for b in range(3)] for a in range(3)]
        Z = self.zeros()
        Rarc, Sarc = Z, Z
        full = R.X
        t = 0
        for _ in range(self.VX + full + 6):
            t += 1
            R.cap = min(t + 1, full)
            Rvalid = Z if zero_m else [[R.sub(Rarc[a][b], self.Rarc0[a][b]) for b in range(3)] for a in range(3)]
            VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
            E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
            newR = self.rarc(E, lam)
            newS = self.sarc(RealNoFace, Rvalid, chord, lam)
            if R.cap >= full:
                same = all(R.equal(newR[a][b], Rarc[a][b]) and R.equal(newS[a][b], Sarc[a][b])
                           for a in range(3) for b in range(3))
            else:
                same = False
            Rarc, Sarc = newR, newS
            if same:
                break
        else:
            raise RuntimeError('arc fixed point did not converge')
        R.cap = full
        if zero_m:
            self.Rarc0 = Rarc
            Rvalid = Z
        else:
            Rvalid = [[R.sub(Rarc[a][b], self.Rarc0[a][b]) for b in range(3)] for a in range(3)]
        VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
        E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
        child = matadd(R, matadd(R, VirtR, Sarc), ParcR)
        return dict(E=E, Real=Real, RealNoFace=RealNoFace, Rvalid=Rvalid, chord=chord, child=child, lam=lam)

    def virtuals(self, Rvalid, Sarc):
        R = self.R
        ok = lambda a, b: not (ISV[a] and ISV[b])
        VirtR = [[R.mul(self.Fin, Rvalid[a][b]) if ok(a, b) else None for b in range(3)] for a in range(3)]
        chord = [[R.add(VirtR[a][b], Sarc[a][b]) for b in range(3)] for a in range(3)]
        ParcR = [[R.mul(self.Fin, chord[a][b]) if ok(a, b) else None for b in range(3)] for a in range(3)]
        return VirtR, chord, ParcR

    def rarc(self, E, lam):
        R = self.R
        LE = [[R.mul(lam[i], E[i][j]) for j in range(3)] for i in range(3)]
        out = self.zeros()
        cur = E
        for d in range(3, self.dmax + 1):
            cur = matmul(R, cur, LE)
            if d in self.cX:
                for a in range(3):
                    for b in range(3):
                        out[a][b] = R.add(out[a][b], R.mul(self.cX[d], cur[a][b]))
            if all(cur[a][b] is None for a in range(3) for b in range(3)):
                break
        return out

    def edge_weight(self, et, l, b, RealNoFace, Rvalid, chord):
        return RealNoFace[l][b] if et == 'real' else (Rvalid[l][b] if et == 'R' else chord[l][b])

    def sarc(self, RealNoFace, Rvalid, chord, lam):
        R = self.R
        out = self.zeros()
        for a in range(3):
            res = self.walk(a, RealNoFace, Rvalid, chord, lam)
            for b in range(3):
                out[a][b] = R.mul(self.Fin, res[b])
        return out

    def walk(self, a, RealNoFace, Rvalid, chord, lam):
        R = self.R
        states = {(a, False, 1 if ISV[a] else 0, 0): Mono(0, 0)}
        result = {b: None for b in range(3)}
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, vs, ne), ser in states.items():
                for et in ('real', 'R', 'P'):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and not ANCHOR[l]:
                        continue
                    for b in range(3):
                        w = self.edge_weight(et, l, b, RealNoFace, Rvalid, chord)
                        if w is None:
                            continue
                        if ISV[b]:
                            if vs == 0:
                                nvs = 1
                            elif vs == 1 and real:
                                nvs = 3
                            else:
                                continue
                        else:
                            nvs = {0: 0, 1: 2, 2: 2, 3: 3}[vs]
                        s2 = R.mul(ser, w)
                        if s2 is None:
                            continue
                        nne = min(ne + 1, 2)
                        if nne >= 2:
                            result[b] = R.add(result[b], s2)
                        key = (b, real, nvs, nne)
                        nxt = R.mul(s2, lam[b])
                        if nxt is not None:
                            new[key] = R.add(new.get(key), nxt)
            states = {k: v for k, v in new.items() if v is not None}
            if not states:
                break
        return result

    # ------------------------------------------------------------------ rooted totals
    def root(self, A):
        R = self.R
        E, Real, child, lam = A['E'], A['Real'], A['child'], A['lam']
        anchors = [0, 1]
        total = None
        RL = [[R.mul(Real[i][j], lam[j]) for j in range(3)] for i in range(3)]
        EL = [[R.mul(E[i][j], lam[j]) for j in range(3)] for i in range(3)]
        cur = matmul(R, RL, EL)
        for d in range(3, self.dmax + 1):
            close = matmul(R, cur, E)
            if d in self.cX:
                for rho in anchors:
                    total = R.add(total, R.mul(R.mul(self.cX[d], lam[rho]), close[rho][rho]))
            cur = matmul(R, cur, EL)
            if all(cur[a][b] is None for a in range(3) for b in range(3)):
                break
        for d in range(3, self.dmax + 1):
            if d not in self.cX:
                continue
            for rho in anchors:
                path = R.mul(lam[rho], Real[rho][P_])
                for _ in range(d - 3):
                    path = R.mul(R.mul(path, lam[P_]), Real[P_][P_])
                path = R.mul(path, lam[P_])
                for u in range(3):
                    s = R.mul(R.mul(R.mul(path, Real[P_][u]), lam[u]), child[u][rho])
                    total = R.sub(total, R.mul(self.cX[d], s))
        for rho in anchors:
            cyc = self.cycle(rho, A)
            total = R.add(total, R.mul(R.mul(self.Fin, lam[rho]), cyc))
        return total

    def cycle(self, rho, A):
        R = self.R
        RealNoFace, Rvalid, chord, lam = A['RealNoFace'], A['Rvalid'], A['chord'], A['lam']
        firstV = bool(ISV[rho])
        states = {(rho, True, 1 if firstV else 0, False, firstV, False, 0): Mono(0, 0)}
        out = None
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, nV, pairok, lastV, pending, ne), ser in states.items():
                for et in (('real',) if ne == 0 else ('real', 'R', 'P')):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and not ANCHOR[l]:
                        continue
                    for b in range(3):
                        w = self.edge_weight(et, l, b, RealNoFace, Rvalid, chord)
                        if w is None:
                            continue
                        s2 = R.mul(ser, w)
                        if s2 is None:
                            continue
                        if b == rho and ne + 1 >= 3:
                            ok = True
                            if pending:
                                ok = lastV and real
                            elif nV == 2 and not pairok:
                                ok = False
                            if ok:
                                out = R.add(out, s2)
                        if pending:
                            continue
                        nnV, npair, npend = nV, pairok, False
                        if ISV[b]:
                            if nV == 0:
                                nnV = 1
                            elif nV == 1 and lastV and real:
                                nnV, npair = 2, True
                            elif nV == 1 and firstV:
                                nnV, npend = 2, True
                            else:
                                continue
                        key = (b, real, nnV, npair, bool(ISV[b]), npend, min(ne + 1, 3))
                        nxt = R.mul(s2, lam[b])
                        if nxt is not None:
                            new[key] = R.add(new.get(key), nxt)
            states = {k: v for k, v in new.items() if v is not None}
            if not states:
                break
        return out


def rooted_at(G, pts):
    """Rooted totals at marker points pts (n x 3 int array), as a host array (n, VX+1, FY+1)."""
    m = tuple(cp.asarray(np.ascontiguousarray(pts[:, k]), dtype=cp.int64) % G.p for k in range(3))
    r = G.root(G.arcs(m))
    return G.R.to_host(r, len(pts))
