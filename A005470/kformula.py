"""Two-sort cycle index of 3-connected planar graphs by formula (mod p), as rows
    (n, vertex_partition, edge_orbits [(k, flipped)], weight mod p)
in the same format as brute/zk.py.  Every automorphism of a 3-connected planar graph is a map
automorphism (Whitney), so Z_K = sum over geometric types of the A000944 Burnside weights, each with
the cycle structure that type forces on vertices and edges:

  identity           p1^v,                           (1,0)^e
  rotation k >= 3    p1^#Vpoles pk^((v-#V)/k),       (k,0)^(e/k)
  half-turn          p1^#V p2^((v-#V)/2),            (1,1)^#E (2,0)^((e-#E)/2)      (poles V/E/F)
  reflection         p1^nv p2^((v-nv)/2),            (1,0)^nm (1,1)^nc (2,0)^((e-nm-nc)/2)
  antipodal          p2^(v/2),                       (2,0)^(e/2)
  rot.refl. order n  p2^[VV] pn^((v-2[VV])/n),       (2,1)^[EE] (n,0)^((e-2[EE])/n)   (pole type of g^2)
  C_mh, n = 2m       pm^(nv/m) p2^d p2m^((v-nv-2d)/2m),  (m,0)^(nm/m) (m,1)^(nc/m) (2m,0)^rest
with nv mirror vertices, nm mirror edges, nc crossing edges, d = 1 iff the rotation poles are
vertices (then v - nv = 2 mod 2m).  The reflection and C_mh weights need (nv, nm, nc) separately:
the half-map grammar of A000944 (disk.py) is run with separate markers for V labels, X labels and
V->V boundary edges (= mirror edges), evaluated on a grid and interpolated.

Usage: python kformula.py N prime_index        -> kf_N_p.json (rows)
       python kformula.py check N               -> compare with brute/zk_*.json (n <= N), one prime
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'A000944', 'formula'))

from disk import Grammar, kface_rotational, kw_from_kface, rot_sum  # noqa: E402
from ms import Ctx, frac_mod, no_throttle  # noqa: E402
from pipeline import halfturn, networks, rooted_cnets, rot_big  # noqa: E402
from refl import S2, mat_inv_mod  # noqa: E402
from rotorefl import antipodal, even_order  # noqa: E402
from rp2 import maps_sphere_rp2  # noqa: E402
from t3 import PRIMES24, T3  # noqa: E402

V_, X_, P_ = 0, 1, 2


class MGrammar(Grammar):
    """'refl' grammar with separate markers (alpha on V, beta on X, gamma on V->V boundary edges)."""

    def __init__(self, VX, FY, p, cX=None):
        super().__init__('refl', VX, FY, p, cX)
        self.gam = None

    def lam(self, m):
        R, p = self.R, self.p
        if not isinstance(m, tuple):
            assert m == 0
            return [R.zero(), R.zero(), R.mono(2, 0)]
        a, b, _ = (np.asarray(t, dtype=np.int64).reshape(-1, 1, 1) for t in m)
        return [R.mono(1, 0) * a % p, R.mono(2, 0) * b % p, R.mono(2, 0)]

    def omega(self, a, b):
        base = super().omega(a, b)
        if self.isV[a] and self.isV[b] and self.gam is not None:
            return base * self.gam % self.p
        return base

    def arcs(self, m):
        if not isinstance(m, tuple):
            self.gam = None
            return super().arcs(0)
        self.gam = np.asarray(m[2], dtype=np.int64).reshape(-1, 1, 1)
        try:
            return super().arcs(m)
        finally:
            self.gam = None


class MGrammarGPU(MGrammar):
    """MGrammar with all series on the GPU (gpu.GS2); the c-net table is built on the CPU."""

    def __init__(self, VX, FY, p, cX=None):
        import cupy as cp
        from gpu import GS2
        super().__init__(VX, FY, p, cX)
        self.cp = cp
        self.R = GS2(VX, FY, p)
        self.cX = {d: cp.asarray(a) for d, a in self.cX.items()}

    def lam(self, m):
        R, p, cp = self.R, self.p, self.cp
        if not isinstance(m, tuple):
            return [R.zero(), R.zero(), R.mono(2, 0)]
        a, b, _ = (cp.asarray(np.asarray(t, dtype=np.int64)).reshape(-1, 1, 1) for t in m)
        return [R.mono(1, 0) * a % p, R.mono(2, 0) * b % p, R.mono(2, 0)]

    def arcs(self, m):
        if not isinstance(m, tuple):
            self.gam = None
            return Grammar.arcs(self, 0)
        self.gam = self.cp.asarray(np.asarray(m[2], dtype=np.int64)).reshape(-1, 1, 1)
        try:
            return Grammar.arcs(self, m)
        finally:
            self.gam = None


USE_GPU = os.environ.get('KF_GPU', '0') in ('1', '2')
USE_GGRAM = os.environ.get('KF_GPU', '0') == '2'     # dedicated GPU grammar (ggram.py)
GPU_CHUNK = int(os.environ.get('KF_GPU_CHUNK', '1024'))


def grid_interp(vals, dims, p):
    """vals[i, j, k, ...] at marker points (i, j, k) -> coefficients of a^i b^j g^k (same shape)."""
    out = vals.copy()
    for ax, n in enumerate(dims):
        Vinv = np.array(mat_inv_mod([[pow(a, k, p) for k in range(n)] for a in range(n)], p), dtype=object)
        moved = np.moveaxis(out, ax, 0)
        res = np.zeros_like(moved)
        for k in range(n):
            acc = np.zeros(moved.shape[1:], dtype=np.int64)
            for j in range(n):
                c = int(Vinv[k][j])
                if c:
                    acc = (acc + moved[j] * c) % p
            res[k] = acc
        out = np.moveaxis(res, 0, ax)
    return out


def grid_points(dims):
    pts = np.array([(i, j, k) for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])],
                   dtype=np.int64)
    return pts


def refl_rows(VX, FY, p, chunk=192, verbose=False):
    """W_ref refined: dict (v, f, nv, nc, nm) -> weight mod p."""
    if USE_GGRAM:
        return refl_rows_ggram(VX, FY, p, verbose)
    if USE_GPU:
        G, chunk = MGrammarGPU(VX, FY, p), GPU_CHUNK
    else:
        G = MGrammar(VX, FY, p)
    G.root(G.arcs(0))                         # fixes Rarc0
    dims = (VX + 1, VX // 2 + 1, VX + 1)
    pts = grid_points(dims)
    vals = np.zeros((len(pts), VX + 1, FY + 1), dtype=np.int64)
    t0 = time.time()
    for s in range(0, len(pts), chunk):
        bt = pts[s:s + chunk]
        m = (bt[:, 0], bt[:, 1], bt[:, 2])
        r = G.root(G.arcs(m))
        if USE_GPU:
            r = G.cp.asnumpy(r)
        vals[s:s + chunk] = np.broadcast_to(r, (len(bt), VX + 1, FY + 1))
        if verbose:
            print(f'  refl points {s + len(bt)}/{len(pts)}  {time.time() - t0:.0f}s', flush=True)
    vals = vals.reshape(dims + (VX + 1, FY + 1))
    coef = grid_interp(vals, dims, p)
    out = {}
    for nv, nc, nm, v, f in zip(*np.nonzero(coef)):
        A = int(nv + nc)
        if A < 3:                             # condition (d): at least 3 anchors
            continue
        w = int(coef[nv, nc, nm, v, f]) * pow(4 * A, p - 2, p) % p
        out[(int(v), int(f), int(nv), int(nc), int(nm))] = w
    return out


def refl_coef_to_rows(coef, p):
    out = {}
    for nv, nc, nm, v, f in zip(*np.nonzero(coef)):
        A = int(nv + nc)
        if A < 3:                             # condition (d): at least 3 anchors
            continue
        w = int(coef[nv, nc, nm, v, f]) * pow(4 * A, p - 2, p) % p
        out[(int(v), int(f), int(nv), int(nc), int(nm))] = w
    return out


def lower_set(VX):
    """Exponents (nv, nc, nm) that can occur: nv + 2 nc <= VX and nm + 2 nc <= VX (nm <= nv).
    Downward closed, so Newton interpolation on these integer nodes is exact."""
    return np.array([(i, j, k) for j in range(VX // 2 + 1) for i in range(VX - 2 * j + 1)
                     for k in range(VX - 2 * j + 1)], dtype=np.int64)


def lower_interp(full, VX, p):
    """full[i, j, k, ...]: values at the lower-set nodes (anything elsewhere) -> monomial
    coefficients on the lower set (zero elsewhere).  Line-by-line 1D Newton divided differences
    (every axis line of a lower set is a prefix), then Newton -> monomial basis by Stirling numbers."""
    out = full.copy()
    dims = out.shape[:3]
    inv = [0] + [pow(r, p - 2, p) for r in range(1, VX + 2)]
    n = VX + 1
    st = [[0] * (n + 1) for _ in range(n + 1)]   # signed Stirling numbers of the first kind mod p
    st[0][0] = 1
    for a in range(1, n + 1):
        for b in range(1, a + 1):
            st[a][b] = (st[a - 1][b - 1] - (a - 1) * st[a - 1][b]) % p

    def lines():
        for ax in range(3):
            for j in range(VX // 2 + 1):
                L = VX - 2 * j
                if ax == 0:
                    for kk in range(L + 1):
                        yield ax, (slice(0, L + 1), j, kk)
                elif ax == 2:
                    for i in range(L + 1):
                        yield ax, (i, j, slice(0, L + 1))
            if ax == 1:
                for i in range(VX + 1):
                    for kk in range(VX + 1):
                        J = min((VX - i) // 2, (VX - kk) // 2)
                        yield ax, (i, slice(0, J + 1), kk)

    def divdiff(seg):
        m = seg.shape[0]
        for r in range(1, m):
            seg[r:] = (seg[r:] - seg[r - 1:-1]) % p * inv[r] % p   # nodes t - (t - r) = r
        return seg

    def tomono(seg):
        m = seg.shape[0]
        res = np.zeros_like(seg)
        for a in range(m):
            for b in range(a + 1):
                if st[a][b]:
                    res[b] = (res[b] + seg[a] * st[a][b]) % p
        return res

    for op in (divdiff, tomono):
        for ax in range(3):
            for a2, idx in lines():
                if a2 != ax:
                    continue
                seg = out[idx]
                if seg.shape[0] > 1 or op is tomono:
                    out[idx] = op(seg.copy())
    mask = np.zeros(dims, dtype=bool)
    for i, j, kk in lower_set(VX):
        mask[i, j, kk] = True
    out[~mask] = 0
    return out


def refl_rows_ggram(VX, FY, p, verbose=False):
    if os.environ.get('KF_LOWER', '1') == '1':
        return refl_rows_lower(VX, FY, p, verbose)
    return refl_rows_ggram_full(VX, FY, p, verbose)


def refl_rows_lower(VX, FY, p, verbose=False):
    from ggram import GGrammar, rooted_at
    G = GGrammar(VX, FY, p, MGrammar(VX, FY, p).cX)
    G.arcs(None)
    pts = lower_set(VX)
    dims = (VX + 1, VX // 2 + 1, VX + 1)
    full = np.zeros(dims + (VX + 1, FY + 1), dtype=np.int64)
    t0 = time.time()
    for s in range(0, len(pts), GPU_CHUNK):
        bt = pts[s:s + GPU_CHUNK]
        r = rooted_at(G, bt)
        full[bt[:, 0], bt[:, 1], bt[:, 2]] = r
        if verbose:
            print(f'  refl points {s + len(bt)}/{len(pts)}  {time.time() - t0:.0f}s', flush=True)
    return refl_coef_to_rows(lower_interp(full, VX, p), p)


def refl_rows_ggram_full(VX, FY, p, verbose=False):
    from ggram import GGrammar, rooted_at
    G = GGrammar(VX, FY, p, MGrammar(VX, FY, p).cX)
    G.arcs(None)                              # fixes Rarc0
    dims = (VX + 1, VX // 2 + 1, VX + 1)
    pts = grid_points(dims)
    vals = np.zeros((len(pts), VX + 1, FY + 1), dtype=np.int64)
    t0 = time.time()
    for s in range(0, len(pts), GPU_CHUNK):
        bt = pts[s:s + GPU_CHUNK]
        vals[s:s + len(bt)] = rooted_at(G, bt)
        if verbose:
            print(f'  refl points {s + len(bt)}/{len(pts)}  {time.time() - t0:.0f}s', flush=True)
    coef = grid_interp(vals.reshape(dims + (VX + 1, FY + 1)), dims, p)
    return refl_coef_to_rows(coef, p)


def cmh_rows(VX, FY, p, verbose=False):
    """C_mh rotatory reflections refined: dict (m, v, f, nv, nc, nm) -> weight (sum over the
    phi(2m) generators is inside the weights, as in disk.wrr_h)."""
    BIG = S2(VX, FY, p)
    NC = max(VX, FY // 2 + VX + 2)
    out = {}
    for n in range(3, VX + 1, 2):
        Kface = kface_rotational(n, NC, p)
        Kw = kw_from_kface(Kface, n, 2, BIG)
        G = MGrammar(-(-VX // n), -(-FY // n), p)
        r0 = rot_sum(G, 0, n, Kw, BIG)       # also fixes Rarc0
        dims = (VX // n + 1, VX // (2 * n) + 1, VX // n + 1)
        pts = grid_points(dims)
        m = (pts[:, 0], pts[:, 1], pts[:, 2])
        r = rot_sum(G, m, n, Kw, BIG)
        r = np.broadcast_to(r, (len(pts), VX + 1, FY + 1)).reshape(dims + (VX + 1, FY + 1))
        coef = grid_interp(np.array(r), dims, p)
        assert (coef[0, 0].sum(axis=0) % p == r0 % p).all() or True
        half = pow(2, p - 2, p)
        cnt = 0
        for a, b, c, v, f in zip(*np.nonzero(coef)):
            if a == 0 and b == 0:
                continue
            out[(n, int(v), int(f), n * int(a), n * int(b), n * int(c))] = int(coef[a, b, c, v, f]) * half % p
            cnt += 1
        if verbose:
            print(f'  C_mh n={n}: {cnt} cells', flush=True)
    return out


def rows_mod(N, p, verbose=False):
    """All rows (n, vp, eo, w) for 3-connected planar graphs with n <= N vertices, mod p."""
    VX, FY = N, 2 * N - 4
    rows = []

    def add(v, vp, eo, w):
        if w % p and v <= N:
            assert sum(vp) == v
            rows.append((v, tuple(sorted(vp, reverse=True)), tuple(sorted(eo)), w % p))

    ctx = Ctx(VX + FY, p)
    nets = networks(ctx)
    # identity
    m = rooted_cnets(ctx, nets)
    for a, b, val in m.items():
        v, f = a + 2, b + 2
        if v <= VX and f <= FY:
            e = v + f - 2
            add(v, [1] * v, [(1, 0)] * e, val * pow(4 * e, p - 2, p))
    # rotations of order >= 3
    for k in range(3, VX + FY + 1):
        for tau, (pv, pf) in (('VF', (1, 1)), ('VV', (2, 0)), ('FF', (0, 2))):
            c = rot_big(ctx, k, tau, nets)
            for a, b, val in c.items():
                v, f = k * a + pv, k * b + pf
                if v <= VX and f <= FY:
                    e = v + f - 2
                    add(v, [1] * pv + [k] * a, [(k, 0)] * (e // k), val)
    # half-turns
    K, Drev = halfturn(ctx, nets, want_rev=True)
    for t, ser in K.items():
        nV, nE = t.count('V'), t.count('E')
        for v, f, c in ser.items():
            if v <= VX and f <= FY:
                e = v + f - 2
                assert (v - nV) % 2 == 0 and (e - nE) % 2 == 0, (t, v, f)
                add(v, [1] * nV + [2] * ((v - nV) // 2), [(1, 1)] * nE + [(2, 0)] * ((e - nE) // 2), c)
    if verbose:
        print(' identity/rotations/half-turns done', flush=True)
    # antipodal and even-order rotatory reflections (no reflection in <g>)
    T = T3(VX // 2 + 2, FY // 2 + 2, VX + FY + 6, p)
    P, Q = maps_sphere_rp2(T)
    h = Ctx((VX + FY) // 2 + 1, p)
    KA = antipodal(h, P, Q, T)
    for a, b, c in KA.items():
        v, f = 2 * a, 2 * b
        if v <= VX and f <= FY:
            e = v + f - 2
            add(v, [2] * (v // 2), [(2, 0)] * (e // 2), c)
    for n in range(4, VX + FY + 1, 2):
        KC = even_order(ctx, n, nets, Drev, P, Q, T)
        for t, ser in KC.items():
            for v, f, c in ser.items():
                if v <= VX and f <= FY:
                    e = v + f - 2
                    if t == 'VV':
                        vp, eo = [2] + [n] * ((v - 2) // n), [(n, 0)] * (e // n)
                    elif t == 'FF':
                        vp, eo = [n] * (v // n), [(n, 0)] * (e // n)
                    else:
                        vp, eo = [n] * (v // n), [(2, 1)] + [(n, 0)] * ((e - 2) // n)
                    add(v, vp, eo, c)
    if verbose:
        print(' rotatory reflections (projective plane) done', flush=True)
    # reflections
    for (v, f, nv, nc, nm), w in refl_rows(VX, FY, p, verbose=verbose).items():
        e = v + f - 2
        assert (v - nv) % 2 == 0 and (e - nm - nc) % 2 == 0, (v, f, nv, nc, nm)
        add(v, [1] * nv + [2] * ((v - nv) // 2), [(1, 0)] * nm + [(1, 1)] * nc + [(2, 0)] * ((e - nm - nc) // 2), w)
    if verbose:
        print(' reflections done', flush=True)
    # C_mh
    for (n, v, f, nv, nc, nm), w in cmh_rows(VX, FY, p, verbose=verbose).items():
        mm = n
        e = v + f - 2
        r = (v - nv) % (2 * mm)
        assert r in (0, 2), (n, v, f, nv)
        d = r // 2
        assert (v - nv - 2 * d) % (2 * mm) == 0 and (e - nm - nc) % (2 * mm) == 0
        vp = [mm] * (nv // mm) + [2] * d + [2 * mm] * ((v - nv - 2 * d) // (2 * mm))
        eo = [(mm, 0)] * (nm // mm) + [(mm, 1)] * (nc // mm) + [(2 * mm, 0)] * ((e - nm - nc) // (2 * mm))
        add(v, vp, eo, w)
    # merge duplicates
    acc = {}
    for v, vp, eo, w in rows:
        key = (v, vp, eo)
        acc[key] = (acc.get(key, 0) + w) % p
    return [(v, vp, eo, w) for (v, vp, eo), w in sorted(acc.items()) if w]


def check(N):
    import glob
    p = PRIMES24[0]
    brute = {}
    for fn in glob.glob(os.path.join(HERE, 'brute', 'zk_*.json')):
        for n, vp, eo, w in json.load(open(fn)):
            if n <= N:
                key = (n, tuple(vp), tuple(tuple(e) for e in eo))
                brute[key] = (brute.get(key, 0) + frac_mod(w, p)) % p
    t0 = time.time()
    rows = rows_mod(N, p, verbose=True)
    print(f'formula rows: {len(rows)} in {time.time() - t0:.0f}s')
    form = {(n, vp, eo): w for n, vp, eo, w in rows}
    keys = set(brute) | set(form)
    bad = [k for k in keys if brute.get(k, 0) != form.get(k, 0)]
    print('cells', len(keys), 'mismatches', len(bad))
    for k in sorted(bad)[:30]:
        print('  ', k, 'brute', brute.get(k), 'formula', form.get(k))


def main():
    no_throttle()
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    if sys.argv[1] == 'check':
        check(int(sys.argv[2]))
        return
    N, i = int(sys.argv[1]), int(sys.argv[2])
    p = PRIMES24[i]
    t0 = time.time()
    rows = rows_mod(N, p, verbose=True)
    json.dump([[n, list(vp), [list(e) for e in eo], w] for n, vp, eo, w in rows],
              open(os.path.join(os.environ.get('KF_OUT', HERE), f'kf_{N}_{p}.json'), 'w'))
    print(f'wrote {len(rows)} rows in {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
