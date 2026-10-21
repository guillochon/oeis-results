"""Generic disk grammar (refl.py generalized) and its rotationally symmetric version.

A "disk" is a 2-connected plane graph H with simple outer cycle C' whose cone is 3-connected, with
labels on C'. Two modes:
  'refl': labels V, X, P (mirror vertex, crossing-edge end, plain): the half of a reflective c-net,
          weights in the doubled map M (x^2 per strict vertex, y^2 per inner face, y per s-piece).
  'cone': labels A, P: A = neighbour of an extra vertex w; H + w is a c-net X (rooted at a vertex):
          weights in X (x per vertex, y per face; each A starts one face at w).
In both modes the 3-connectivity conditions are the same local rules (refl.py docstring), except that
the "V rule" (at most one mirror vertex per face, mirror edges between consecutive V) exists only in
'refl' mode.

rot_sum(): Burnside weight of disks with a rotation rho of odd order n >= 3 about an interior point:
    sum over (H, labels, rho) of 1/|centralizer| = (fixed centre node) x (one period of the boundary).
The centre is an S-node (a face touching C' in n-periodic points) or an R-node: a c-net with a
rotation of order n whose pole is its outer face, counted by outer degree (Kface). Kface comes from
the 'cone' mode by duality (a rotation pole face of X is a pole vertex of its dual), recursively.
"""
from math import gcd

import numpy as np

from refl import S2, mat_inv_mod, matadd, matmul
from t3 import T3, cnets_by_rootface


def phi(n):
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


class Grammar:
    def __init__(self, mode, VX, FY, p, cX=None):
        self.mode, self.VX, self.FY, self.p = mode, VX, FY, p
        self.R = R = S2(VX, FY, p)
        if mode == 'refl':
            self.anchor, self.isV, self.lexp, self.iv = [1, 1, 0], [1, 0, 0], [1, 2, 2], 2
        else:
            self.anchor, self.isV, self.lexp, self.iv = [1, 0], [0, 0], [1, 1], 1
        self.L = len(self.anchor)
        self.plain = self.L - 1
        self.Fin = R.mono(0, self.iv)
        iv = self.iv
        if cX is None:
            T = T3(VX, FY // iv + VX + 1, VX, p)
            c, _ = cnets_by_rootface(T)
            cX = {}
            for v in range(4, T.X + 1):
                for f in range(4, T.Y + 1):
                    for d in range(3, min(v, T.Z) + 1):
                        val = int(c[v, f, d])
                        i, j = iv * (v - d), iv * (f - 1 - d)
                        if val and i <= VX and j <= FY:
                            cX.setdefault(d, R.zero())
                            cX[d][i, j] = (cX[d][i, j] + val) % p
        self.cX = cX
        self.dmax = max(cX) if cX else 3
        self.Rarc0 = None

    def zeros(self):
        return [[self.R.zero() for _ in range(self.L)] for _ in range(self.L)]

    def omega(self, a, b):
        R = self.R
        if self.anchor[a] and not (self.isV[a] and self.isV[b]):
            return R.mono(0, 1)
        return R.mono(0, 0)

    def pair_ok(self, a, b):
        return not (self.isV[a] and self.isV[b])

    def lam(self, alpha):
        """alpha: an int, or an array of values evaluated as one batch (leading axis)."""
        R, p = self.R, self.p
        if isinstance(alpha, np.ndarray):
            alpha = alpha.reshape(-1, 1, 1)
        return [R.mono(e, 0) * (alpha if self.anchor[i] else 1) % p for i, e in enumerate(self.lexp)]

    # ------------------------------------------------------------------ arcs
    def arcs(self, alpha):
        """Arc fixed point at anchor marker alpha; alpha = 0 must be computed first (stores Rarc0).
        alpha may be an array of nonzero values: one batched evaluation (leading axis)."""
        R, p, L = self.R, self.p, self.L
        zero_alpha = not isinstance(alpha, np.ndarray) and alpha == 0
        lam = self.lam(alpha)
        Real = [[R.mul(self.omega(a, b), self.Fin) for b in range(L)] for a in range(L)]
        RealNoFace = [[self.omega(a, b) for b in range(L)] for a in range(L)]
        Z = self.zeros()
        Rarc, Sarc = Z, Z
        for _ in range(self.VX + 4):
            Rvalid = Z if zero_alpha else [[(Rarc[a][b] - self.Rarc0[a][b]) % p for b in range(L)] for a in range(L)]
            VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
            E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
            newR = self.rarc(E, lam)
            newS = self.sarc(RealNoFace, Rvalid, chord, lam)
            same = all((newR[a][b] == Rarc[a][b]).all() and (newS[a][b] == Sarc[a][b]).all()
                       for a in range(L) for b in range(L))
            Rarc, Sarc = newR, newS
            if same:
                break
        else:
            raise RuntimeError('arc fixed point did not converge')
        if zero_alpha:
            self.Rarc0 = Rarc
            Rvalid = Z
        else:
            Rvalid = [[(Rarc[a][b] - self.Rarc0[a][b]) % p for b in range(L)] for a in range(L)]
        VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
        E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
        child = matadd(R, matadd(R, VirtR, Sarc), ParcR)
        return dict(E=E, Real=Real, RealNoFace=RealNoFace, Rvalid=Rvalid, chord=chord, child=child, lam=lam)

    def virtuals(self, Rvalid, Sarc):
        R, p, L = self.R, self.p, self.L
        VirtR = [[R.mul(self.Fin, Rvalid[a][b]) if self.pair_ok(a, b) else R.zero() for b in range(L)]
                 for a in range(L)]
        chord = [[(VirtR[a][b] + Sarc[a][b]) % p for b in range(L)] for a in range(L)]
        ParcR = [[R.mul(self.Fin, chord[a][b]) if self.pair_ok(a, b) else R.zero() for b in range(L)]
                 for a in range(L)]
        return VirtR, chord, ParcR

    def rarc(self, E, lam):
        R, p, L = self.R, self.p, self.L
        LE = [[R.mul(lam[i], E[i][j]) for j in range(L)] for i in range(L)]
        out = self.zeros()
        cur = E
        for d in range(3, self.dmax + 1):
            cur = matmul(R, cur, LE)
            if d in self.cX:
                for a in range(L):
                    for b in range(L):
                        out[a][b] = (out[a][b] + R.mul(self.cX[d], cur[a][b])) % p
            if not any(cur[a][b].any() for a in range(L) for b in range(L)):
                break
        return out

    def edge_weight(self, et, l, b, RealNoFace, Rvalid, chord):
        return RealNoFace[l][b] if et == 'real' else (Rvalid[l][b] if et == 'R' else chord[l][b])

    def sarc(self, RealNoFace, Rvalid, chord, lam):
        R, p, L = self.R, self.p, self.L
        out = self.zeros()
        for a in range(L):
            res = self.walk(a, RealNoFace, Rvalid, chord, lam)
            for b in range(L):
                out[a][b] = R.mul(self.Fin, res[b])
        return out

    def walk(self, a, RealNoFace, Rvalid, chord, lam):
        R, p, L = self.R, self.p, self.L
        states = {(a, False, 1 if self.isV[a] else 0, 0): R.mono(0, 0)}
        result = {b: R.zero() for b in range(L)}
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, vs, ne), ser in states.items():
                for et in ('real', 'R', 'P'):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and not self.anchor[l]:
                        continue
                    for b in range(L):
                        w = self.edge_weight(et, l, b, RealNoFace, Rvalid, chord)
                        if not w.any():
                            continue
                        if self.isV[b]:
                            if vs == 0:
                                nvs = 1
                            elif vs == 1 and real:
                                nvs = 3
                            else:
                                continue
                        else:
                            nvs = {0: 0, 1: 2, 2: 2, 3: 3}[vs]
                        s2 = R.mul(ser, w)
                        nne = min(ne + 1, 2)
                        if nne >= 2:
                            result[b] = (result[b] + s2) % p
                        key = (b, real, nvs, nne)
                        new[key] = (new.get(key, R.zero()) + R.mul(s2, lam[b])) % p
            states = new
            if not states:
                break
        return result

    # ------------------------------------------------------------------ rooted totals (no symmetry)
    def root(self, A):
        R, p, L = self.R, self.p, self.L
        E, Real, child, lam = A['E'], A['Real'], A['child'], A['lam']
        anchors = [i for i in range(L) if self.anchor[i]]
        P_ = self.plain
        total = R.zero()
        RL = [[R.mul(Real[i][j], lam[j]) for j in range(L)] for i in range(L)]
        EL = [[R.mul(E[i][j], lam[j]) for j in range(L)] for i in range(L)]
        cur = matmul(R, RL, EL)
        for d in range(3, self.dmax + 1):
            close = matmul(R, cur, E)
            if d in self.cX:
                for rho in anchors:
                    total = (total + R.mul(R.mul(self.cX[d], lam[rho]), close[rho][rho])) % p
            cur = matmul(R, cur, EL)
            if not any(cur[a][b].any() for a in range(L) for b in range(L)):
                break
        for d in range(3, self.dmax + 1):
            if d not in self.cX:
                continue
            for rho in anchors:
                path = R.mul(lam[rho], Real[rho][P_])
                for _ in range(d - 3):
                    path = R.mul(R.mul(path, lam[P_]), Real[P_][P_])
                path = R.mul(path, lam[P_])
                for u in range(L):
                    s = R.mul(R.mul(R.mul(path, Real[P_][u]), lam[u]), child[u][rho])
                    total = (total - R.mul(self.cX[d], s)) % p
        for rho in anchors:
            cyc = self.cycle(rho, A)
            total = (total + R.mul(R.mul(self.Fin, lam[rho]), cyc)) % p
        return total

    def cycle(self, rho, A):
        R, p, L = self.R, self.p, self.L
        RealNoFace, Rvalid, chord, lam = A['RealNoFace'], A['Rvalid'], A['chord'], A['lam']
        firstV = bool(self.isV[rho])
        states = {(rho, True, 1 if firstV else 0, False, firstV, False, 0): R.mono(0, 0)}
        out = R.zero()
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, nV, pairok, lastV, pending, ne), ser in states.items():
                for et in (('real',) if ne == 0 else ('real', 'R', 'P')):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and not self.anchor[l]:
                        continue
                    for b in range(L):
                        w = self.edge_weight(et, l, b, RealNoFace, Rvalid, chord)
                        if not w.any():
                            continue
                        s2 = R.mul(ser, w)
                        if b == rho and ne + 1 >= 3:
                            ok = True
                            if pending:
                                ok = lastV and real
                            elif nV == 2 and not pairok:
                                ok = False
                            if ok:
                                out = (out + s2) % p
                        if pending:
                            continue
                        nnV, npair, npend = nV, pairok, False
                        if self.isV[b]:
                            if nV == 0:
                                nnV = 1
                            elif nV == 1 and lastV and real:
                                nnV, npair = 2, True
                            elif nV == 1 and firstV:
                                nnV, npend = 2, True
                            else:
                                continue
                        key = (b, real, nnV, npair, bool(self.isV[b]), npend, min(ne + 1, 3))
                        new[key] = (new.get(key, R.zero()) + R.mul(s2, lam[b])) % p
            states = new
            if not states:
                break
        return out


def interpolate(vals, alphas, p, zero):
    n = len(alphas)
    Vinv = mat_inv_mod([[pow(a, k, p) for k in range(n)] for a in alphas], p)
    out = []
    for k in range(n):
        acc = zero.copy()
        for j in range(n):
            if Vinv[k][j]:
                acc = (acc + vals[j] * int(Vinv[k][j])) % p
        out.append(acc)
    return out


def rooted_by_anchors(G):
    """[alpha^A] of the anchor-rooted totals, A = 0..VX."""
    alphas = list(range(G.VX + 1))
    vals = [G.root(G.arcs(a)) for a in alphas]
    return interpolate(vals, alphas, G.p, G.R.zero())


# ---------------------------------------------------------------------- rotational centre
def scale_up(a, n, BIG):
    """small array a[..., i, j] -> big array with entry at (..., n i, n j)."""
    out = np.zeros(a.shape[:-2] + (BIG.VX + 1, BIG.FY + 1), dtype=np.int64)
    r, c = min(a.shape[-2], BIG.VX // n + 1), min(a.shape[-1], BIG.FY // n + 1)
    out[..., :n * (r - 1) + 1:n, :n * (c - 1) + 1:n] = a[..., :r, :c]
    return out


def rot_parts(G, beta, n, BIG):
    """The Kw-independent parts of rot_sum: ({d': trace of (E Lam)^d' per period}, S-centre total)."""
    p, L = G.p, G.L
    A = G.arcs(beta)
    up = lambda M: [[scale_up(M[a][b], n, BIG) for b in range(L)] for a in range(L)]
    E, RNF, Rv, ch, lam = up(A['E']), up(A['RealNoFace']), up(A['Rvalid']), up(A['chord']), A['lam']
    lamb = [scale_up(l, n, BIG) for l in lam]
    # R-centre: trace of (E Lam)^d' per period
    trs = {}
    EL = [[BIG.mul(E[i][j], lamb[j]) for j in range(L)] for i in range(L)]
    cur = EL
    for dp in range(1, BIG.VX + 1):
        tr = BIG.zero()
        for a in range(L):
            tr = (tr + cur[a][a]) % p
        if tr.any():
            trs[dp] = tr
        cur = matmul(BIG, cur, EL)
        if not any(cur[a][b].any() for a in range(L) for b in range(L)):
            break
    # S-centre: periodic chain, no mirror vertex among the touches; state (label, edge into it real)
    total = BIG.zero()
    Fin = BIG.mono(0, G.iv)          # the centre face is a single face (not n copies)
    states = [(l, r) for l in range(L) if not G.isV[l] for r in (True, False)]
    idx = {s: i for i, s in enumerate(states)}
    T = [[BIG.zero() for _ in states] for _ in states]
    for (l, r) in states:
        for et in ('real', 'R', 'P'):
            real = et == 'real'
            if r and real and not G.anchor[l]:
                continue
            for b in range(L):
                if G.isV[b]:
                    continue
                w = RNF[l][b] if real else (Rv[l][b] if et == 'R' else ch[l][b])
                if w.any():
                    j = idx[(b, real)]
                    T[idx[(l, r)]][j] = (T[idx[(l, r)]][j] + BIG.mul(w, lamb[b])) % p
    cur = T
    for kp in range(1, BIG.VX + 1):
        if n * kp >= 3:
            tr = BIG.zero()
            for i in range(len(states)):
                tr = (tr + cur[i][i]) % p
            if tr.any():
                coef = phi(n) * pow(2 * n * kp, p - 2, p) % p
                total = (total + BIG.mul(Fin, tr) * coef) % p
        cur = matmul(BIG, cur, T)
        if not any(c.any() for row in cur for c in row):
            break
    return trs, total


def rot_sum(G, beta, n, Kw, BIG, parts=None):
    """sum over (disk, labels, rotation of order n) of 1/|centralizer|, as a BIG series.
    G: grammar truncated at (BIG.VX // n, BIG.FY // n); beta: marker per anchor of one period.
    Kw[d']: BIG series of R-centre weights for outer degree n d' (sum of Kface times inner weights).
    parts: rot_parts(G, beta, n, BIG), if already computed (it does not depend on Kw)."""
    trs, total = parts or rot_parts(G, beta, n, BIG)
    for dp, tr in trs.items():
        if dp in Kw:
            total = (total + BIG.mul(Kw[dp], tr)) % G.p
    return total


def memo_arcs(G):
    cache = {}
    orig = G.arcs

    def arcs(alpha):
        key = tuple(alpha.tolist()) if isinstance(alpha, np.ndarray) else alpha
        if key not in cache:
            if 0 not in cache and key != 0:
                cache[0] = orig(0)
            cache[key] = orig(alpha)
        return cache[key]
    G.arcs = arcs
    return G


def kw_from_kface(Kface, n, iv, BIG):
    """R-centre weights: Kw[d'] = sum_{v,f} Kface[n d'][v, f] x^(iv (v - n d')) y^(iv (f - 1 - n d'))."""
    Kw = {}
    for d, arr in Kface.items():
        if d % n:
            continue
        acc = BIG.zero()
        for v, f in zip(*np.nonzero(arr)):
            i, j = iv * (v - d), iv * (f - 1 - d)
            if 0 <= i <= BIG.VX and 0 <= j <= BIG.FY:
                acc[i, j] = (acc[i, j] + arr[v, f]) % BIG.p
        if acc.any():
            Kw[d // n] = acc
    return Kw


def kface_rotational(n, NC, p, verbose=False):
    """Kface[d][v, f]: sum over c-nets X (v vertices, f faces) with a designated face O of degree d and
    a rotation of order n (odd >= 3) about O, of 1/|Aut(X, O)|.  Square truncation v, f <= NC."""
    BIG = S2(NC, NC, p)
    G = memo_arcs(Grammar('cone', -(-NC // n), -(-NC // n), p))
    amax = NC // n + 1
    betas = list(range(amax + 1))
    p0 = rot_parts(G, 0, n, BIG)                                    # beta = 0 first (fixes Rarc0)
    pb = rot_parts(G, np.array(betas[1:], dtype=np.int64), n, BIG)  # the others as one batch
    Kface = {}
    for it in range(NC + 2):
        Kw = kw_from_kface(Kface, n, 1, BIG)
        vb = np.broadcast_to(rot_sum(G, None, n, Kw, BIG, pb), (len(betas) - 1, BIG.VX + 1, BIG.FY + 1))
        vals = [rot_sum(G, 0, n, Kw, BIG, p0)] + [vb[i] for i in range(len(betas) - 1)]
        coef = interpolate(vals, betas, p, BIG.zero())
        new = {}
        for a in range(1, len(coef)):
            if coef[a].any():
                kv = np.zeros_like(coef[a])
                kv[1:, :] = coef[a][:-1, :]          # the cone vertex w
                new[n * a] = kv.T.copy()             # dual: pole vertex of degree d -> pole face
        same = set(new) == set(Kface) and all(np.array_equal(new[d], Kface[d]) for d in new)
        Kface = new
        if verbose:
            print('iteration', it, sorted(Kface), flush=True)
        if same:
            return Kface
    raise RuntimeError('Kface recursion did not converge')


def wrr_h(VX, FY, p, nmax=None, verbose=False):
    """W_rr from C_nh symmetry (mirror x rotation of odd order n >= 3), as an S2 array [v, f] in M:
    1/2 * (rotationally symmetric labelled halves with at least one anchor per period)."""
    BIG = S2(VX, FY, p)
    total = BIG.zero()
    NC = max(VX, FY // 2 + VX + 2)
    for n in range(3, (nmax or VX) + 1, 2):
        Kface = kface_rotational(n, NC, p)
        Kw = kw_from_kface(Kface, n, 2, BIG)
        G = memo_arcs(Grammar('refl', -(-VX // n), -(-FY // n), p))
        r1 = rot_sum(G, 1, n, Kw, BIG)
        r0 = rot_sum(G, 0, n, Kw, BIG)
        part = (r1 - r0) % p * pow(2, p - 2, p) % p
        if verbose:
            print('n', n, 'nonzero cells', int(np.count_nonzero(part)), flush=True)
        total = (total + part) % p
    return total
