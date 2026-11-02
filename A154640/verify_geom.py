"""Independent check of all seven sequences by plain floating-point geometry.

Nothing here uses the integer sphere vectors, the reflection matrices R_i or the Coxeter group.
Objects are real spheres (center, radius) and planes (unit normal n, offset h: n.x = h).
Mallows's start: three unit spheres centered at (2/sqrt3, 0, 0), (-1/sqrt3, +-1, 0); planes z = +-1.

Descartes reflection of member E of a quintuple {A, B, C, D, E}: invert E in the sphere through the
tangency points of A, B, C, D (it is orthogonal to all four, so the inversion fixes them and maps E to
the other sphere tangent to all four).  Spheres are identified by rounding (bend, center).

Strategy (a): BFS over ordered configurations; a configuration's generation is its BFS depth
              (= word length); a sphere's generation is the first depth at which it appears.
Strategy (b): a quintuple's generation is that of its youngest sphere; generation g+1 spheres are the
              new spheres made by reflecting generation-g quintuples (closing level g first).
Every quintuple produced is checked to consist of five mutually tangent objects.
"""
import math
import sys
from collections import Counter

import numpy as np

S3 = math.sqrt(3)
TOL = 1e-7


class Sph:
    __slots__ = ("c", "r")
    def __init__(self, c, r): self.c = np.asarray(c, float); self.r = float(r)


class Pl:
    __slots__ = ("n", "h")
    def __init__(self, n, h): self.n = np.asarray(n, float); self.h = float(h)


def bend(o): return 0.0 if isinstance(o, Pl) else 1.0 / o.r


def key(o):
    if isinstance(o, Pl):
        n, h = (o.n, o.h) if o.h >= 0 else (-o.n, -o.h)
        return ("P",) + tuple(round(v * 1e6) for v in n) + (round(h * 1e6),)
    b = 1.0 / o.r
    bi = round(b)
    assert abs(b - bi) < 1e-6 * max(1, b), ("bend not integral", b)
    return (bi,) + tuple(round(v * 1e6) for v in o.c)


def tangent_point(p, q):
    if isinstance(p, Pl) and isinstance(q, Pl):
        return None  # parallel planes touch at infinity
    if isinstance(p, Pl):
        p, q = q, p
    if isinstance(q, Pl):
        return p.c - (q.n @ p.c - q.h) * q.n * 1.0  # foot of the perpendicular (|n.c - h| = r)
    d = q.c - p.c
    return p.c + d * (p.r / np.linalg.norm(d))


def tangent(p, q):
    if isinstance(p, Pl) and isinstance(q, Pl):
        return abs(abs(p.n @ q.n) - 1) < 1e-9
    if isinstance(p, Pl):
        p, q = q, p
    if isinstance(q, Pl):
        return abs(abs(q.n @ p.c - q.h) - p.r) < TOL * max(1, p.r)
    return abs(np.linalg.norm(p.c - q.c) - (p.r + q.r)) < TOL * max(1, p.r + q.r)


def dual_sphere(four):
    """Sphere through the tangency points of four mutually tangent objects, or
    ("plane", n, h) when the two plates are among them (the points are then coplanar)."""
    pts = []
    for i in range(4):
        for j in range(i + 1, 4):
            t = tangent_point(four[i], four[j])
            if t is not None:
                pts.append(t)
    P = np.array(pts)
    _, sv, vt = np.linalg.svd(P - P.mean(axis=0))
    if sv[-1] < 1e-9 * max(1, sv[0]):
        n = vt[-1] / np.linalg.norm(vt[-1])
        return "plane", n, float(n @ P.mean(axis=0))
    # |x|^2 + a.x + c = 0 through the points: least squares on all (5 or 6) points
    A = np.array([[*p, 1.0] for p in pts])
    rhs = -np.array([p @ p for p in pts])
    sol, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    a, c = sol[:3], sol[3]
    O = -a / 2
    rho2 = O @ O - c
    resid = max(abs(np.linalg.norm(p - O) ** 2 - rho2) for p in pts)
    assert rho2 > 0 and resid < 1e-7 * max(1, rho2), ("dual sphere fit", resid, rho2)
    return O, rho2


def invert(o, O, rho2):
    if isinstance(o, Pl):
        d = o.h - o.n @ O                       # signed distance from O to the plane
        r = rho2 / (2 * abs(d))
        return Sph(O + o.n * (rho2 / (2 * d)), r)
    v = o.c - O
    den = v @ v - o.r ** 2
    if abs(den) < 1e-9 * max(1, o.r ** 2):     # sphere through O: image is a plane
        u = v / np.linalg.norm(v)
        return Pl(u, u @ O + rho2 / (2 * o.r))
    k = rho2 / den
    return Sph(O + k * v, abs(k) * o.r)


def reflect(objs, i):
    four = [objs[k] for k in range(5) if k != i]
    d = dual_sphere(four)
    if isinstance(d[0], str):
        _, n, h = d
        o = objs[i]
        assert isinstance(o, Sph)
        return Sph(o.c - 2 * (n @ o.c - h) * n, o.r)   # mirror image
    return invert(objs[i], *d)


START = [Sph((2 / S3, 0, 0), 1), Sph((-1 / S3, 1, 0), 1), Sph((-1 / S3, -1, 0), 1),
         Pl((0, 0, 1), 1), Pl((0, 0, -1), 1)]


class Store:
    def __init__(self):
        self.obj, self.ids, self.gen = [], {}, []

    def add(self, o, g):
        k = key(o)
        if k in self.ids:
            return self.ids[k], False
        self.ids[k] = len(self.obj); self.obj.append(o); self.gen.append(g)
        return self.ids[k], True


def check_quint(store, q):
    objs = [store.obj[i] for i in q]
    for i in range(5):
        for j in range(i + 1, 5):
            assert tangent(objs[i], objs[j]), "reflected quintuple is not mutually tangent"


def strategy_a(K):
    st = Store()
    q0 = tuple(st.add(o, 0)[0] for o in START)
    seen = {q0}; frontier = [q0]
    out = [Counter()]
    for i in q0:
        out[0]["n"] += 1; out[0]["s"] += bend(st.obj[i]); out[0]["t"] += bend(st.obj[i]) ** 2
    for k in range(1, K + 1):
        nxt = []
        c = Counter()
        for q in frontier:
            objs = [st.obj[i] for i in q]
            for i in range(5):
                new = reflect(objs, i)
                sid, isnew = st.add(new, k)
                q2 = q[:i] + (sid,) + q[i + 1:]
                if q2 in seen:
                    continue
                seen.add(q2); nxt.append(q2)
                if isnew:
                    b = bend(new); c["n"] += 1; c["s"] += b; c["t"] += b * b
        out.append(c)
        frontier = nxt
        print(f"  (a) generation {k}: {len(nxt)} configurations, {c['n']} new spheres", file=sys.stderr)
    return out, st


def strategy_b(K, check=True):
    st = Store()
    q0 = tuple(st.add(o, 0)[0] for o in START)
    seen = {frozenset(q0)}
    level = [q0]
    out = [Counter({"n": 5, "s": 3.0, "t": 3.0})]
    quints = []
    for g in range(K):
        stack = list(level); nq = len(level); nxt = []
        c = Counter()
        while stack:
            q = stack.pop()
            objs = [st.obj[i] for i in q]
            for i in range(5):
                new = reflect(objs, i)
                sid, isnew = st.add(new, g + 1)
                if isnew:
                    b = bend(new); c["n"] += 1; c["s"] += b; c["t"] += b * b
                q2 = q[:i] + (sid,) + q[i + 1:]
                k2 = frozenset(q2)
                if k2 in seen:
                    continue
                seen.add(k2)
                if check:
                    check_quint(st, q2)
                gq = max(st.gen[s] for s in q2)
                assert gq >= g, "quintuple found after its own generation was closed"
                if gq == g:
                    stack.append(q2); nq += 1
                else:
                    nxt.append(q2)
        quints.append(nq)
        out.append(c)
        level = nxt
        print(f"  (b) generation {g} closed: {nq} quintuples; generation {g + 1}: {c['n']} spheres", file=sys.stderr)
    return out, quints, st


if __name__ == "__main__":
    Ka = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    Kb = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    oa, _ = strategy_a(Ka)
    print("A154640", [int(c["n"]) for c in oa])
    print("A154641", [round(c["s"] / 3, 3) for c in oa])
    print("A154642", [round(c["t"] / 3, 3) for c in oa])
    ob, quints, _ = strategy_b(Kb)
    print("quintuples(b)", quints)
    print("A154643", [int(c["n"]) for c in ob])
    print("A154644", [round(c["s"] / 3, 3) for c in ob])
    print("A154645", [round(c["t"] / 3, 3) for c in ob])
