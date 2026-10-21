"""Reflection-symmetric c-nets: W_ref(v,f) = sum over (M, reflection) of 1/|centralizer|.

A reflective c-net M <-> its half H (2-connected plane graph, simple outer cycle C', vertices of C'
labelled V = mirror vertex, X = end of a crossing edge, P = plain), M = double of H. M is
3-connected iff (a) the cone over C' is a c-net, (b) every inner face of H has <= 1 V or exactly two
V joined by a mirror edge on it, (c) every inner face meets every piece of C' (arc between
consecutive anchors V/X) in nothing, a vertex or an edge; plus >= 3 anchors (brute-checked in
brute/refl_check.py).  W_ref = 1/2 * sum_{(H,lab)} 1/|Aut(H,lab)|.

Counting: H decomposes along its 2-cuts (all on C') into R-nodes (c-nets whose outer face carries all
virtual edges) and S-nodes (polygons = faces of H touching C' >= 3 times), no P-nodes, no S-S
adjacency. A virtual edge of an R-node carries an R-arc (2-touch face) or an S-arc; S-nodes carry
R-arcs. Conditions become local (each pocket = R-arc has an interior anchor; touches of an S-face
between two real chain edges are anchors; V-count per multi-touch face) except one correction at the
root.  Rooted at an anchor rho and the C'-edge leaving it; anchors marked by alpha, evaluated at
alpha = 0..VX and interpolated; W_ref = 1/2 sum_A [alpha^A] rooted / (2A), A >= 3.
Labels index 0 = V, 1 = X, 2 = P.
"""
import sys

import numpy as np

from fastmul import mulmod
from t3 import PRIMES24, T3, cnets_by_rootface

V_, X_, P_ = 0, 1, 2


class S2:
    """2D series mod p with rectangular truncation (x <= VX, y <= FY)."""

    def __init__(self, VX, FY, p):
        self.VX, self.FY, self.p = VX, FY, p

    def zero(self):
        return np.zeros((self.VX + 1, self.FY + 1), dtype=np.int64)

    def mono(self, i, j, c=1):
        a = self.zero()
        if i <= self.VX and j <= self.FY:
            a[i, j] = c % self.p
        return a

    def mul(self, a, b):
        return mulmod(a, b, self.p, (self.VX + 1, self.FY + 1))


def matmul(R, A, B):
    n, m, q = len(A), len(B), len(B[0])
    out = [[R.zero() for _ in range(q)] for _ in range(n)]
    for i in range(n):
        for k in range(m):
            if not A[i][k].any():
                continue
            for j in range(q):
                if B[k][j].any():
                    out[i][j] = (out[i][j] + R.mul(A[i][k], B[k][j])) % R.p
    return out


def matadd(R, A, B):
    return [[(a + b) % R.p for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def matscale_cols(R, A, lam):
    """A * diag(lam)."""
    return [[R.mul(A[i][j], lam[j]) for j in range(3)] for i in range(len(A))]


class Refl:
    def __init__(self, VX, FY, p):
        self.R = R = S2(VX, FY, p)
        self.VX, self.FY, self.p = VX, FY, p
        # c-nets by root face degree, mapped to H-weights x^(2(v-d)) y^(2(f-1-d))
        T = T3(VX, FY // 2 + VX + 1, VX, p)
        c, _ = cnets_by_rootface(T)
        self.cX = {}
        for v in range(4, T.X + 1):
            for f in range(4, T.Y + 1):
                for d in range(3, min(v, T.Z) + 1):
                    val = int(c[v, f, d])
                    if not val:
                        continue
                    i, j = 2 * (v - d), 2 * (f - 1 - d)
                    if i <= VX and j <= FY:
                        self.cX.setdefault(d, R.zero())
                        self.cX[d][i, j] = (self.cX[d][i, j] + val) % p
        self.dmax = max(self.cX) if self.cX else 3

    def omega(self, a, b):
        """C'-edge a -> b: y for the first edge of an s-piece (leaves an anchor, not V->V)."""
        R = self.R
        if a in (V_, X_) and not (a == V_ and b == V_):
            return R.mono(0, 1)
        return R.mono(0, 0)

    def run(self, alpha):
        """alpha: 0 (first; sets Rarc0), or an array of nonzero values evaluated as one batch."""
        R, p = self.R, self.p
        x, x2, y2 = R.mono(1, 0), R.mono(2, 0), R.mono(0, 2)
        batched = isinstance(alpha, np.ndarray)
        if batched:
            alpha = alpha.reshape(-1, 1, 1)
        else:
            assert alpha == 0
        lam = [x * alpha % p, x2 * alpha % p, x2]
        Real = [[R.mul(self.omega(a, b), y2) for b in range(3)] for a in range(3)]
        RealNoFace = [[self.omega(a, b) for b in range(3)] for a in range(3)]
        Z = [[R.zero() for _ in range(3)] for _ in range(3)]
        Rarc = Z
        Sarc = Z
        Rarc0 = self.Rarc0 if batched else None
        for it in range(self.VX + 3):
            if batched:
                Rvalid = [[(Rarc[a][b] - Rarc0[a][b]) % p for b in range(3)] for a in range(3)]
            else:
                Rvalid = Z  # anchor-free R-arcs never satisfy the pocket condition
            VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
            E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
            newR = self.rarc(E, lam)
            newS = self.sarc(RealNoFace, Rvalid, chord, lam)
            same = all((newR[a][b] == Rarc[a][b]).all() and (newS[a][b] == Sarc[a][b]).all()
                       for a in range(3) for b in range(3))
            Rarc, Sarc = newR, newS
            if same:
                break
        else:
            raise RuntimeError('arc fixed point did not converge')
        if not batched:
            self.Rarc0 = Rarc
            Rvalid = Z
        else:
            Rvalid = [[(Rarc[a][b] - Rarc0[a][b]) % p for b in range(3)] for a in range(3)]
        VirtR, chord, ParcR = self.virtuals(Rvalid, Sarc)
        E = matadd(R, matadd(R, matadd(R, Real, VirtR), Sarc), ParcR)
        child = matadd(R, matadd(R, VirtR, Sarc), ParcR)
        return self.root(E, Real, child, RealNoFace, Rvalid, chord, lam)

    def virtuals(self, Rvalid, Sarc):
        """Pocket options on a virtual edge a -> b.
        VirtR: R-arc under an R-node (2-touch face spanning both, y^2, not both ends V).
        chord: P-node = chord a-b + child arc (R-arc: closed 2-touch face y^2, not both V; or S-arc,
               which then needs an interior anchor - automatic, S-arcs without anchors vanish).
        ParcR: chord arc under an R-node: the R-node's face closed by the chord is a 2-touch face."""
        R, p = self.R, self.p
        y2 = R.mono(0, 2)
        nvv = lambda a, b: not (a == V_ and b == V_)
        VirtR = [[R.mul(y2, Rvalid[a][b]) if nvv(a, b) else R.zero() for b in range(3)] for a in range(3)]
        chord = [[(VirtR[a][b] + Sarc[a][b]) % p for b in range(3)] for a in range(3)]
        ParcR = [[R.mul(y2, chord[a][b]) if nvv(a, b) else R.zero() for b in range(3)] for a in range(3)]
        return VirtR, chord, ParcR

    def rarc(self, E, lam):
        """sum_d cX_d * E (Lam E)^(d-2)."""
        R, p = self.R, self.p
        LE = [[R.mul(lam[i], E[i][j]) for j in range(3)] for i in range(3)]
        out = [[R.zero() for _ in range(3)] for _ in range(3)]
        cur = E  # E (Lam E)^(k)
        for d in range(3, self.dmax + 1):
            cur = matmul(R, cur, LE)  # E (LE)^(d-2)
            if d in self.cX:
                for a in range(3):
                    for b in range(3):
                        out[a][b] = (out[a][b] + R.mul(self.cX[d], cur[a][b])) % p
            if not any(cur[a][b].any() for a in range(3) for b in range(3)):
                break
        return out

    def sarc(self, RealNoFace, Rvalid, chord, lam):
        """S-arc a -> b: y^2 * sum over chains a = z0, ..., zk = b, k >= 2 (closed by the parent's
        virtual edge, which can never be part of the allowed V pair)."""
        R, p = self.R, self.p
        y2 = R.mono(0, 2)
        out = [[R.zero() for _ in range(3)] for _ in range(3)]
        for a in range(3):
            res = self.walk(a, RealNoFace, Rvalid, chord, lam)
            for b in range(3):
                out[a][b] = R.mul(y2, res[b])
        return out

    def walk(self, a, RealNoFace, Rvalid, chord, lam):
        """Chains from label a; returns b -> series over chains with >= 2 edges ending at label b
        (end vertex unweighted). (c2): an intermediate vertex between two real edges is an anchor.
        V rule: <= 1 V among all chain vertices, or exactly 2 adjacent via a real chain edge.
        vs: 0 = no V, 1 = one V and it is the last vertex, 2 = one V not last, 3 = pair done."""
        R, p = self.R, self.p
        states = {(a, False, 1 if a == V_ else 0, 0): R.mono(0, 0)}
        result = {b: R.zero() for b in range(3)}
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, vs, ne), ser in states.items():
                for et in ('real', 'R', 'P'):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and l == P_:
                        continue
                    for b in range(3):
                        w = RealNoFace[l][b] if real else (Rvalid[l][b] if et == 'R' else chord[l][b])
                        if not w.any():
                            continue
                        if b == V_:
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

    def root(self, E, Real, child, RealNoFace, Rvalid, chord, lam):
        R, p = self.R, self.p
        total = R.zero()
        LE = [[R.mul(lam[i], E[i][j]) for j in range(3)] for i in range(3)]
        # Root-R: sum_d cX_d sum_{rho anchor} lam_rho [Real Lam (E Lam)^(d-2) E]_{rho,rho}
        RL = [[R.mul(Real[i][j], lam[j]) for j in range(3)] for i in range(3)]
        EL = [[R.mul(E[i][j], lam[j]) for j in range(3)] for i in range(3)]
        cur = matmul(R, RL, EL)  # Real Lam (E Lam)^(d-2) for d = 3
        for d in range(3, self.dmax + 1):
            # closing with * E -> d edges in total
            close = matmul(R, cur, E)
            if d in self.cX:
                for rho in (V_, X_):
                    total = (total + R.mul(R.mul(self.cX[d], lam[rho]), close[rho][rho])) % p
            cur = matmul(R, cur, EL)
            if not any(cur[a][b].any() for a in range(3) for b in range(3)):
                break
        # correction: last edge (u_{d-1}, rho) a child containing all other anchors
        for d in range(3, self.dmax + 1):
            if d not in self.cX:
                continue
            for rho in (V_, X_):
                # rho -> P (real) -> P ... -> u_{d-1} (any) : d-1 real edges, d-2 P vertices
                path = R.mul(lam[rho], Real[rho][P_])
                last = P_
                for _ in range(d - 3):
                    path = R.mul(R.mul(path, lam[P_]), Real[P_][P_])
                path = R.mul(path, lam[P_])
                for u in range(3):
                    s = R.mul(R.mul(path, Real[P_][u]), lam[u])
                    s = R.mul(s, child[u][rho])
                    total = (total - R.mul(self.cX[d], s)) % p
        # Root-S: polygon through rho, first edge real; S-face y^2
        y2 = R.mono(0, 2)
        for rho in (V_, X_):
            cyc = self.cycle(rho, RealNoFace, Rvalid, chord, lam)
            total = (total + R.mul(R.mul(y2, lam[rho]), cyc)) % p
        return total

    def cycle(self, rho, RealNoFace, Rvalid, chord, lam):
        """Root polygons (>= 3 edges) through anchor rho, first edge real; (c2) and the V rule
        cyclically. State: (label, lastreal, nV, pairok, lastV, pending, nedges)."""
        R, p = self.R, self.p
        firstV = rho == V_
        states = {(rho, True, 1 if firstV else 0, False, firstV, False, 0): R.mono(0, 0)}
        out = R.zero()
        for _ in range(self.VX + 2):
            new = {}
            for (l, lastreal, nV, pairok, lastV, pending, ne), ser in states.items():
                for et in (('real',) if ne == 0 else ('real', 'R', 'P')):
                    real = et == 'real'
                    if ne >= 1 and lastreal and real and l == P_:
                        continue
                    for b in range(3):
                        w = RealNoFace[l][b] if real else (Rvalid[l][b] if et == 'R' else chord[l][b])
                        if not w.any():
                            continue
                        s2 = R.mul(ser, w)
                        # close at rho (b is rho itself, not a new vertex)
                        if b == rho and ne + 1 >= 3:
                            ok = True
                            if pending:
                                ok = lastV and real
                            elif nV == 2 and not pairok:
                                ok = False
                            if ok:
                                out = (out + s2) % p
                        if pending:
                            continue  # a pending V must be the last vertex
                        # b as a new intermediate vertex
                        nnV, npair, npend = nV, pairok, False
                        if b == V_:
                            if nV == 0:
                                nnV = 1
                            elif nV == 1 and lastV and real:
                                nnV, npair = 2, True
                            elif nV == 1 and firstV:
                                nnV, npend = 2, True   # may pair with rho through the closing edge
                            else:
                                continue
                        key = (b, real, nnV, npair, b == V_, npend, min(ne + 1, 3))
                        new[key] = (new.get(key, R.zero()) + R.mul(s2, lam[b])) % p
            states = new
            if not states:
                break
        return out


def wref(VX, FY, p):
    """W_ref mod p as an array [v, f] (M vertices, faces)."""
    Rf = Refl(VX, FY, p)
    R = Rf.R
    alphas = list(range(0, VX + 1))
    vals = [Rf.run(0)]                    # alpha = 0 first: it fixes the anchor-free arcs Rarc0
    batch = Rf.run(np.array(alphas[1:], dtype=np.int64))
    vals += [batch[i] for i in range(len(alphas) - 1)]
    # interpolate in alpha: values at alpha = 0..VX -> coefficients of alpha^A
    n = len(alphas)
    Vm = np.array([[pow(a, k, p) for k in range(n)] for a in alphas], dtype=object)
    Vinv = mat_inv_mod(Vm, p)
    coeffs = []
    for k in range(n):
        acc = R.zero()
        for j in range(n):
            if Vinv[k][j]:
                acc = (acc + vals[j] * int(Vinv[k][j])) % p
        coeffs.append(acc)
    W = R.zero()
    for A in range(3, n):
        W = (W + coeffs[A] * pow(4 * A, p - 2, p)) % p
    return W


def mat_inv_mod(M, p):
    n = len(M)
    A = [[int(M[i][j]) % p for j in range(n)] + [int(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        r = next(r for r in range(c, n) if A[r][c])
        A[c], A[r] = A[r], A[c]
        iv = pow(A[c][c], p - 2, p)
        A[c] = [v * iv % p for v in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [(v - f * w) % p for v, w in zip(A[r], A[c])]
    return [row[n:] for row in A]


if __name__ == '__main__':
    import glob
    import json
    from fractions import Fraction
    from ms import frac_mod
    VX, FY = int(sys.argv[1]), int(sys.argv[2])
    p = PRIMES24[0]
    W = wref(VX, FY, p)
    stats = {}
    for fn in glob.glob('../brute/stats_*.json'):
        for key, d in json.load(open(fn)).items():
            stats.setdefault(key, {}).update(d)
    bad = ok = 0
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v > VX or f > FY:
            continue
        want = sum((Fraction(w) for t, w in d.items() if json.loads(t)[0] == 'ref'), Fraction(0))
        if int(W[v, f]) == frac_mod(want, p):
            ok += 1
        else:
            bad += 1
            print('W_ref mismatch', v, f, 'brute', want)
    print('ok', ok, 'bad', bad)
