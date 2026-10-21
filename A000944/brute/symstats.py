"""Brute-force symmetry statistics of polyhedra (3-connected planar graphs = c-nets).

Reads graph6 files of 2-connected planar graphs (gen_poly.sh), keeps the 3-connected ones,
builds the (unique up to mirror) embedding, enumerates all map automorphisms and classifies
each non-identity one by geometric type:
  ('rot', k, poles)        rotation of order k, poles = sorted pair of cell types V/E/F
  ('ref', nv, nm, nc)      reflection: nv mirror vertices, nm mirror edges, nc crossing edges
  ('rr', order)            rotatory reflection (no fixed cell; order 2 = antipodal map)
Output: JSON {"v,f": {type: weight}} with weight = sum over maps of (#elements of type)/|Aut|,
so the per-(v,f) sum over all types (identity included) is the number of polyhedra U(v,f).

Usage: python symstats.py planar2c_*.g6 > stats.json
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction

import networkx as nx


def embed(G):
    ok, emb = nx.check_planarity(G)
    assert ok
    darts = [(u, w) for u in G for w in emb.neighbors_cw_order(u)]
    idx = {d: i for i, d in enumerate(darts)}
    n = len(darts)
    sig = [0] * n
    for u in G:
        nb = list(emb.neighbors_cw_order(u))
        for i, w in enumerate(nb):
            sig[idx[(u, w)]] = idx[(u, nb[(i + 1) % len(nb)])]
    alp = [idx[(w, u)] for (u, w) in darts]
    return darts, sig, alp


def automorphisms(sig, alp):
    n = len(sig)
    sinv = [0] * n
    for i, j in enumerate(sig):
        sinv[j] = i
    out = []
    for rev in (False, True):
        s2 = sinv if rev else sig
        for d in range(n):
            m = [-1] * n
            m[0] = d
            stack = [0]
            ok = True
            while stack and ok:
                x = stack.pop()
                for a, b in ((sig[x], s2[m[x]]), (alp[x], alp[m[x]])):
                    if m[a] < 0:
                        m[a] = b
                        stack.append(a)
                    elif m[a] != b:
                        ok = False
                        break
            if ok and len(set(m)) == n:
                out.append((rev, m))
    return out


def cells(darts, sig, alp):
    n = len(darts)
    phi = [sig[alp[i]] for i in range(n)]
    seen = [False] * n
    faces = []
    for i in range(n):
        if not seen[i]:
            f = []
            j = i
            while not seen[j]:
                seen[j] = True
                f.append(j)
                j = phi[j]
            faces.append(f)
    return faces


def classify(rev, m, darts, alp, faces, vert_of, face_of):
    n = len(m)
    order = 1
    p = m[:]
    while p != list(range(n)):
        p = [m[x] for x in p]
        order += 1
    fixv = {vert_of[i] for i in range(n) if vert_of[m[i]] == vert_of[i]}
    fixf = {face_of[i] for i in range(n) if face_of[m[i]] == face_of[i]} if not rev else None
    mirror = sum(1 for i in range(n) if m[i] == i) // 2
    cross = sum(1 for i in range(n) if m[i] == alp[i]) // 2
    if not rev:
        if rev is False and fixf is not None:
            # orientation-preserving: faces map as dart cycles
            pass
        poles = ['V'] * len(fixv) + ['E'] * cross + ['F'] * len(fixf)
        assert len(poles) == 2, (order, poles)
        return ('rot', order, ''.join(sorted(poles)))
    if order == 2 and (fixv or mirror or cross):
        return ('ref', len(fixv), mirror, cross)
    assert not fixv and not mirror and not cross
    return ('rr', order)


def face_index(faces, darts, alp, rev_faces=False):
    pass


def main():
    stats = defaultdict(lambda: defaultdict(Fraction))
    for fn in sys.argv[1:]:
        with open(fn, 'rb') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                G = nx.from_graph6_bytes(line)
                if nx.node_connectivity(G) < 3:
                    continue
                darts, sig, alp = embed(G)
                faces = cells(darts, sig, alp)
                vert_of = [u for (u, w) in darts]
                face_of = [0] * len(darts)
                for k, f in enumerate(faces):
                    for d in f:
                        face_of[d] = k
                auts = automorphisms(sig, alp)
                A = len(auts)
                key = f"{G.number_of_nodes()},{len(faces)}"
                for rev, m in auts:
                    if not rev and m == list(range(len(m))):
                        t = ('id',)
                    elif rev:
                        # a reversing map sends face cycles to reversed cycles: compare edge sets
                        t = classify_rev(m, darts, alp, faces, vert_of)
                    else:
                        t = classify(rev, m, darts, alp, faces, vert_of, face_of)
                    stats[key][json.dumps(t)] += Fraction(1, A)
    out = {k: {t: str(w) for t, w in v.items()} for k, v in sorted(stats.items())}
    json.dump(out, sys.stdout, indent=0)


def classify_rev(m, darts, alp, faces, vert_of):
    n = len(m)
    order = 1
    p = m[:]
    while p != list(range(n)):
        p = [m[x] for x in p]
        order += 1
    fixv = {vert_of[i] for i in range(n) if vert_of[m[i]] == vert_of[i]}
    mirror = sum(1 for i in range(n) if m[i] == i) // 2
    cross = sum(1 for i in range(n) if m[i] == alp[i]) // 2
    if order == 2 and (fixv or mirror or cross):
        return ('ref', len(fixv), mirror, cross)
    assert not fixv and not mirror and not cross, (order, fixv, mirror, cross)
    # refine: does <g> contain a reflection (C_nh type)?  else the pole type of g^2
    p = list(range(n))
    powers = []
    for _ in range(order):
        p = [m[x] for x in p]
        powers.append(p)
    half = powers[order // 2 - 1]          # g^(order/2)
    if order % 4 == 2 and order > 2:
        # g^(order/2) has odd exponent: reflection (fixes something) or antipodal map
        if any(vert_of[half[i]] == vert_of[i] for i in range(n)) or any(half[i] in (i, alp[i]) for i in range(n)):
            return ('rr', order, 'h')
    if order == 2:
        return ('rr', 2)
    g2 = powers[1]
    nv = len({vert_of[i] for i in range(n) if vert_of[g2[i]] == vert_of[i]})
    ne = sum(1 for i in range(n) if g2[i] == alp[i]) // 2
    pol = 'V' * (nv // 1) + 'E' * ne
    pol = pol if len(pol) == 2 else (pol + 'F' * (2 - len(pol)))
    return ('rr', order, ''.join(sorted(pol)))


if __name__ == '__main__':
    main()
