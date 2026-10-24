"""Two-disk circle puzzle GG(n, r): pieces cut by the closure of the disk boundaries.

Disks A, B of radius r centred at -1 and +1; each turns by 2*pi/n.  The cut set is the
smallest set of circular arcs that contains both boundary circles and is mapped into itself
by each disk's rotation (applied to the part of the cut set inside that disk).  Every cut is
an arc of a radius-r circle, so we store arcs as (centre, union of angle intervals) and
propagate with a worklist.  Pieces are the bounded faces of the resulting arc arrangement,
counted with Euler's formula F = E - V + C.

    python twodisk.py 5 1.3          # piece count for n = 5, r = 1.3
"""
import cmath
import math
import sys

import numpy as np

TAU = 2 * math.pi
EPS_ANG = 1e-13      # intervals shorter than this are dropped as round-off
EPS_PT = 1e-10        # points / centres closer than this are identified


def norm(t):
    t %= TAU
    return t


def split(s, e):
    """Interval [s, e] (e - s <= TAU) as a list of intervals inside [0, TAU]."""
    if e - s >= TAU - EPS_ANG:
        return [(0.0, TAU)]
    s0 = norm(s)
    e0 = s0 + (e - s)
    if e0 <= TAU:
        return [(s0, e0)]
    return [(s0, TAU), (0.0, e0 - TAU)]


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + EPS_ANG:
            if e > out[-1][1]:
                out[-1] = (out[-1][0], e)
        else:
            out.append((s, e))
    return out


def subtract(iv, have):
    """Parts of the intervals iv not covered by the sorted disjoint list have."""
    out = []
    for s, e in iv:
        cur = s
        for hs, he in have:
            if he <= cur:
                continue
            if hs >= e:
                break
            if hs > cur:
                out.append((cur, min(hs, e)))
            cur = max(cur, he)
            if cur >= e:
                break
        if cur < e:
            out.append((cur, e))
    return [(s, e) for s, e in out if e - s > EPS_ANG]


def intersect(iv, s, e):
    """Intersect sorted intervals iv (in [0, TAU]) with the arc [s, e] (any real s)."""
    out = []
    for a, b in split(s, e):
        for x, y in iv:
            lo, hi = max(a, x), min(b, y)
            if hi - lo > EPS_ANG:
                out.append((lo, hi))
    return out


class Puzzle:
    def __init__(self, n, r, cap=20000):
        self.n, self.r, self.cap = n, r, cap
        self.w = cmath.exp(2j * math.pi / n)
        self.disks = [complex(-1, 0), complex(1, 0)]
        self.centres = []
        self.arcs = []
        self.grid = {}

    def cid(self, c):
        kx, ky = round(c.real / EPS_PT / 10), round(c.imag / EPS_PT / 10)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in self.grid.get((kx + dx, ky + dy), ()):
                    if abs(self.centres[i] - c) < EPS_PT:
                        return i
        i = len(self.centres)
        self.centres.append(c)
        self.arcs.append([])
        self.grid.setdefault((kx, ky), []).append(i)
        return i

    def add(self, i, iv):
        new = subtract(merge(iv), self.arcs[i])
        if new:
            self.arcs[i] = merge(self.arcs[i] + new)
        return new

    def close(self):
        work = []
        for d in self.disks:
            i = self.cid(d)
            work.append((i, self.add(i, [(0.0, TAU)])))
        r = self.r
        step = TAU / self.n
        while work:
            i, iv = work.pop()
            c = self.centres[i]
            for dc in self.disks:
                dist = abs(dc - c)
                if dist < EPS_PT or dist >= 2 * r:
                    continue
                half = math.acos(dist / (2 * r))
                mid = cmath.phase(dc - c)
                inside = intersect(iv, mid - half, mid + half)
                if not inside:
                    continue
                j = self.cid(dc + (c - dc) * self.w)
                img = []
                for s, e in inside:
                    img += split(s + step, e + step)
                new = self.add(j, img)
                if new:
                    work.append((j, new))
                    if len(self.centres) > self.cap:
                        return False
        return True

    def pieces(self):
        r = self.r
        C = np.array(self.centres)
        N = len(C)
        pts_on = [[] for _ in range(N)]   # angles of vertices on circle i

        def on_arcs(i, t):
            t = norm(t)
            for s, e in self.arcs[i]:
                if s - 1e-9 <= t <= e + 1e-9:
                    return True
            return False

        for i in range(N):
            for s, e in self.arcs[i]:
                if e - s < TAU - EPS_ANG:
                    pts_on[i] += [s, e]
            d = C[i + 1:] - C[i]
            dist = np.abs(d)
            for k in np.nonzero(dist < 2 * r - 1e-12)[0]:
                j = i + 1 + k
                h = math.acos(dist[k] / (2 * r))
                m = cmath.phase(d[k])
                for sgn in (1, -1):
                    ti = m + sgn * h
                    p = C[i] + r * cmath.exp(1j * ti)
                    tj = cmath.phase(p - C[j])
                    if on_arcs(i, ti) and on_arcs(j, tj):
                        pts_on[i].append(norm(ti))
                        pts_on[j].append(norm(tj))
        # identify vertices
        vid = {}
        verts = []

        def vkey(p):
            kx, ky = round(p.real / EPS_PT / 10), round(p.imag / EPS_PT / 10)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for v in vid.get((kx + dx, ky + dy), ()):
                        if abs(verts[v] - p) < EPS_PT:
                            return v
            v = len(verts)
            verts.append(p)
            vid.setdefault((kx, ky), []).append(v)
            return v

        edges = []   # (u, v, circle, t0, t1): ccw arc of circle from vertex u at t0 to v at t1
        for i in range(N):
            angs = sorted(set(pts_on[i]))
            for s, e in self.arcs[i]:
                full = e - s >= TAU - EPS_ANG
                inner = [t for t in angs if s - 1e-9 <= t <= e + 1e-9]
                seq = []
                for t in inner:
                    v = vkey(C[i] + r * cmath.exp(1j * t))
                    if not seq or seq[-1][0] != v:
                        seq.append((v, t))
                if full:
                    if len(seq) > 1 and seq[0][0] == seq[-1][0]:
                        seq.pop()
                    if not seq:
                        raise ValueError("isolated circle")
                    seq.append((seq[0][0], seq[0][1] + TAU))
                edges += [(a[0], b[0], i, a[1], b[1]) for a, b in zip(seq, seq[1:])]
        self.edges, self.verts = edges, verts
        V = len(verts)
        parent = list(range(V))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b, *_ in edges:
            parent[find(a)] = find(b)
        comps = len({find(v) for v in range(V)})
        return len(edges) - V + comps

    def pieces_fast(self, delta=0.05, keep_edges=False):
        """pieces() vectorised: candidate arc pairs from a k-d tree over short sub-arcs, then numpy throughout."""
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        from scipy.spatial import cKDTree
        r = self.r
        C = np.array(self.centres)
        ai, ast, aen = [], [], []
        for i, iv in enumerate(self.arcs):
            for s, e in iv:
                ai.append(i); ast.append(s); aen.append(e)
        ai, ast, aen = np.array(ai), np.array(ast), np.array(aen)
        full = aen - ast >= TAU - EPS_ANG
        # sub-arcs of span <= delta, each with a bounding disc (arc midpoint, half the arc length)
        k = np.maximum(1, np.ceil((aen - ast) / delta)).astype(int)
        sub_arc = np.repeat(np.arange(len(ai)), k)
        j = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)
        span = ((aen - ast) / k)[sub_arc]
        s0 = ast[sub_arc] + j * span
        cen = C[ai[sub_arc]] + r * np.exp(1j * (s0 + span / 2))
        rad = r * span / 2 * 1.0001 + 1e-9
        tree = cKDTree(np.c_[cen.real, cen.imag])
        pr = tree.query_pairs(2 * rad.max(), output_type="ndarray")
        a, b = pr[:, 0], pr[:, 1]
        keep = (np.abs(cen[a] - cen[b]) <= rad[a] + rad[b]) & (ai[sub_arc[a]] != ai[sub_arc[b]])
        a, b = a[keep], b[keep]
        ci, cj = C[ai[sub_arc[a]]], C[ai[sub_arc[b]]]
        d = cj - ci
        dist = np.abs(d)
        ok = dist < 2 * r - 1e-12
        a, b, ci, cj, d, dist = a[ok], b[ok], ci[ok], cj[ok], d[ok], dist[ok]
        h = np.arccos(dist / (2 * r))
        m = np.angle(d)
        tol = 1e-9

        def within(t, sub):
            lo = s0[sub] - tol
            return np.mod(t - lo, TAU) <= span[sub] + 2 * tol

        pts_arc, pts_t = [], []
        for sgn in (1, -1):
            ti = m + sgn * h
            p = ci + r * np.exp(1j * ti)
            tj = np.angle(p - cj)
            hit = within(ti, a) & within(tj, b)
            pts_arc += [sub_arc[a[hit]], sub_arc[b[hit]]]
            pts_t += [np.mod(ti[hit], TAU), np.mod(tj[hit], TAU)]
        # clamp angles into their arc's interval (they may sit 1e-9 outside, or wrap past TAU)
        pa = np.concatenate(pts_arc + [np.nonzero(~full)[0]] * 2)
        pt = np.concatenate(pts_t + [ast[~full], aen[~full]])
        wrap = (pt < ast[pa] - 1e-6) & (pt + TAU <= aen[pa] + 1e-6)
        pt = np.where(wrap, pt + TAU, pt)
        pt = np.clip(pt, ast[pa], aen[pa])
        xy = C[ai[pa]] + r * np.exp(1j * pt)
        # identify coincident points
        P = np.c_[xy.real, xy.imag]
        close = cKDTree(P).query_pairs(EPS_PT, output_type="ndarray")
        n_pts = len(P)
        g = coo_matrix((np.ones(len(close)), (close[:, 0], close[:, 1])), shape=(n_pts, n_pts))
        _, vid = connected_components(g, directed=False)
        V = vid.max() + 1
        # walk each arc in angle order
        order = np.lexsort((pt, pa))
        pa, vid, pt = pa[order], vid[order], pt[order]
        dup = np.r_[False, (pa[1:] == pa[:-1]) & (vid[1:] == vid[:-1])]
        pa, vid, pt = pa[~dup], vid[~dup], pt[~dup]
        first = np.r_[True, pa[1:] != pa[:-1]]
        last = np.r_[pa[1:] != pa[:-1], True]
        # closing duplicate on full circles
        fi = np.nonzero(first)[0]; la = np.nonzero(last)[0]
        drop = full[pa[la]] & (vid[la] == vid[fi]) & (la > fi)
        mask = np.ones(len(pa), bool); mask[la[drop]] = False
        pa, vid, pt = pa[mask], vid[mask], pt[mask]
        same = pa[1:] == pa[:-1]
        eu, ev = list(vid[:-1][same]), list(vid[1:][same])
        first = np.r_[True, pa[1:] != pa[:-1]]
        last = np.r_[pa[1:] != pa[:-1], True]
        fi, la = np.nonzero(first)[0], np.nonzero(last)[0]
        fc = full[pa[fi]]
        eu += list(vid[la[fc]]); ev += list(vid[fi[fc]])
        E = len(eu)
        if keep_edges:   # same layout as pieces(): (u, v, circle, t0, t1), for faces()
            circ = ai[pa]
            self.edges = list(zip(eu, ev,
                                  list(circ[:-1][same]) + list(circ[la[fc]]),
                                  list(pt[:-1][same]) + list(pt[la[fc]]),
                                  list(pt[1:][same]) + list(pt[fi[fc]] + TAU)))
        g = coo_matrix((np.ones(E), (eu, ev)), shape=(V, V))
        comps, _ = connected_components(g, directed=False)
        return E - V + comps

    def faces(self):
        """Trace faces of the arrangement; return the areas of the bounded ones."""
        r, C = self.r, self.centres
        out = {}   # vertex -> list of (key, half-edge)
        for k, (u, v, i, t0, t1) in enumerate(self.edges):
            out.setdefault(u, []).append(((norm(t0 + math.pi / 2), 1), 2 * k))
            out.setdefault(v, []).append(((norm(t1 - math.pi / 2), -1), 2 * k + 1))
        pos = {}
        for v, lst in out.items():
            lst.sort()
            for idx, (_, h) in enumerate(lst):
                pos[h] = (v, idx)

        def head(h):
            u, v = self.edges[h >> 1][:2]
            return v if h % 2 == 0 else u

        def area(h):
            _, _, i, t0, t1 = self.edges[h >> 1]
            c = C[i]
            a = 0.5 * (r * c.real * (math.sin(t1) - math.sin(t0))
                       - r * c.imag * (math.cos(t1) - math.cos(t0)) + r * r * (t1 - t0))
            return a if h % 2 == 0 else -a

        seen = [False] * (2 * len(self.edges))
        areas = []
        for h0 in range(len(seen)):
            if seen[h0]:
                continue
            h, a = h0, 0.0
            while not seen[h]:
                seen[h] = True
                a += area(h)
                v, idx = pos[h ^ 1]
                lst = out[v]
                h = lst[(idx - 1) % len(lst)][1]
            areas.append(a)
        return [a for a in areas if a > 0], [a for a in areas if a <= 0]

    def check(self):
        P = self.pieces_fast(keep_edges=True)
        pos, neg = self.faces()
        r = self.r
        union = 2 * math.pi * r * r - (2 * r * r * math.acos(1 / r) - 2 * math.sqrt(r * r - 1))
        return P, len(pos), len(neg), sum(pos) - union, min(pos)


class ExactPuzzle(Puzzle):
    """Puzzle with no accumulated round-off, for large n.

    Centres are kept exactly in Z[zeta_n] (integer vectors modulo the cyclotomic polynomial), so
    they are identified by key rather than by tolerance. Arc endpoints are kept symbolically: a
    label ('c', p, s) is the crossing with circle p at angle arg(C_p - C_i) + s*acos(|C_p - C_i|/2r),
    and ('f', k) is the fixed angle k*2pi/n (the 0 / 2pi seam). A rotation maps labels to labels and
    every angle is recomputed from exact centres, so nothing drifts however long the orbit.
    """
    EPS = 1e-12

    def __init__(self, n, r, cap=20000):
        from sympy import Poly, cyclotomic_poly, symbols
        super().__init__(n, r, cap)
        x = symbols("x")
        phi = Poly(cyclotomic_poly(n, x), x).all_coeffs()[::-1]   # low to high, monic
        self.deg = len(phi) - 1
        self.phi = [int(c) for c in phi[:-1]]
        self.zr = [math.cos(TAU * k / n) for k in range(self.deg)]
        self.zi = [math.sin(TAU * k / n) for k in range(self.deg)]
        self.key = {}
        self.exact = []
        self.iv = []
        one = (1,) + (0,) * (self.deg - 1)
        self.dkeys = [tuple(-c for c in one), one]

    def to_complex(self, v):
        return complex(math.fsum(c * z for c, z in zip(v, self.zr)), math.fsum(c * z for c, z in zip(v, self.zi)))

    def rot(self, v, a):
        """a + (v - a) * zeta, all in Z[zeta]."""
        w = [x - y for x, y in zip(v, a)]
        top = w[-1]
        w = [0] + w[:-1]
        if top:
            w = [x - top * p for x, p in zip(w, self.phi)]
        return tuple(x + y for x, y in zip(w, a))

    def cid(self, v):
        i = self.key.get(v)
        if i is None:
            i = len(self.centres)
            self.key[v] = i
            self.exact.append(v)
            self.centres.append(self.to_complex(v))
            self.iv.append([])
        return i

    def angle(self, i, lab, end):
        if lab[0] == "f":
            t = (lab[1] % self.n) * TAU / self.n
            return TAU if end and t == 0 else t
        d = self.centres[lab[1]] - self.centres[i]
        t = (cmath.phase(d) + lab[2] * math.acos(abs(d) / (2 * self.r))) % TAU
        return TAU if end and t == 0 else t

    def place(self, i, ls, le):
        """Interval(s) on circle i from start label ls to end label le, split at the seam."""
        s, e = self.angle(i, ls, False), self.angle(i, le, True)
        if e >= s:
            return [(s, e, ls, le)]
        return [(s, TAU, ls, ("f", self.n)), (0.0, e, ("f", 0), le)]

    def merge(self, iv):
        out = []
        for it in sorted(iv, key=lambda t: (t[0], t[1])):
            if out and it[0] <= out[-1][1] + self.EPS:
                if it[1] > out[-1][1]:
                    out[-1] = (out[-1][0], it[1], out[-1][2], it[3])
            else:
                out.append(it)
        return out

    def subtract(self, iv, have):
        out = []
        for s, e, ls, le in iv:
            cur, cl = s, ls
            for hs, he, hls, hle in have:
                if he <= cur:
                    continue
                if hs >= e:
                    break
                if hs > cur:
                    out.append((cur, hs, cl, hls) if hs < e else (cur, e, cl, le))
                if he > cur:
                    cur, cl = he, hle
                if cur >= e:
                    break
            if cur < e:
                out.append((cur, e, cl, le))
        return [t for t in out if t[1] - t[0] > self.EPS]

    def add(self, i, iv):
        new = self.subtract(self.merge(iv), self.iv[i])
        if new:
            self.iv[i] = self.merge(self.iv[i] + new)
        return new

    def close(self):
        n, r = self.n, self.r
        work = []
        for v in self.dkeys:
            i = self.cid(v)
            work.append((i, self.add(i, [(0.0, TAU, ("f", 0), ("f", n))])))
        while work:
            i, iv = work.pop()
            for d, dv in enumerate(self.dkeys):
                dc = self.cid(dv)
                if dc == i:
                    continue
                dist = abs(self.centres[dc] - self.centres[i])
                if dist >= 2 * r:
                    continue
                disk = self.place(i, ("c", dc, -1), ("c", dc, 1))
                inside = []
                for a, b, la, lb in disk:
                    for s, e, ls, le in iv:
                        lo, llo = (s, ls) if s >= a else (a, la)
                        hi, lhi = (e, le) if e <= b else (b, lb)
                        if hi - lo > self.EPS:
                            inside.append((llo, lhi))
                if not inside:
                    continue
                j = self.cid(self.rot(self.exact[i], dv))

                def move(lab):
                    if lab[0] == "f":
                        return ("f", lab[1] + 1)
                    return ("c", self.cid(self.rot(self.exact[lab[1]], dv)), lab[2])

                img = []
                for ls, le in inside:
                    img += self.place(j, move(ls), move(le))
                new = self.add(j, img)
                if new:
                    work.append((j, new))
                    if len(self.centres) > self.cap:
                        return False
        self.arcs = [[(s, e) for s, e, _, _ in iv] for iv in self.iv]
        return True


def count(n, r, cap=20000):
    p = Puzzle(n, r, cap)
    if not p.close():
        return None, len(p.centres)
    return p.pieces_fast(), len(p.centres)


if __name__ == "__main__":
    n, r = int(sys.argv[1]), float(sys.argv[2])
    print(count(n, r))
