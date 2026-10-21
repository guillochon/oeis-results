"""Brute-force check of the half-map characterization of reflection-symmetric c-nets.

For every polyhedron W (cone) and vertex w of W whose incident faces are all triangles, H = W - w is a
2-connected plane graph whose outer cycle C' is the link of w.  For every labelling of C' by
V (mirror vertex) / X (crossing-edge end) / P (plain), double H across the mirror:
  V vertices are shared, other vertices get two copies, C'-edges between two V are mirror edges
  (single), every X vertex z gets the crossing edge z_A z_B.
Check: the double M is simple and 3-connected  <=>  conditions (b) and (c) on the inner faces of H.
Also accumulate the Burnside weight sum 1/(2|Stab|) per (v_M, f_M, mirror stats) and compare with the
reflection weights of brute/stats_*.json.

Usage: python refl_check.py NMAX_CONE
"""
import glob
import itertools
import json
import sys
from collections import defaultdict
from fractions import Fraction

import networkx as nx

from symstats import automorphisms, cells, embed


def inner_faces_of_H(darts, faces, w):
    """Faces of W not containing w, as vertex lists (cyclic), = inner faces of H."""
    out = []
    for f in faces:
        vs = [darts[d][0] for d in f]
        if w not in vs:
            out.append(vs)
    return out


def link_cycle(darts, sig, w):
    """Neighbours of w in rotation order."""
    start = next(i for i, d in enumerate(darts) if d[0] == w)
    out = []
    d = start
    while True:
        out.append(darts[d][1])
        d = sig[d]
        if d == start:
            return out


def conditions(lab, link, faces):
    """(b) and (c) for the labelled outer cycle `link` with labels lab[i]."""
    L = len(link)
    pos = {u: i for i, u in enumerate(link)}
    anchors = [i for i in range(L) if lab[i] in 'VX']
    if len(anchors) < 3:   # the mirror must meet M in >= 3 points
        return False
    # piece id of each C'-edge (i, i+1); mirror edges get None
    piece_of_edge = {}
    pid = -1
    start = anchors[0]
    for k in range(L):
        i = (start + k) % L
        if lab[i] in 'VX':
            pid += 1
        j = (i + 1) % L
        if lab[i] == 'V' and lab[j] == 'V':
            piece_of_edge[i] = None
        else:
            piece_of_edge[i] = pid
    # vertex -> set of pieces containing it (anchor belongs to the pieces on both sides)
    pieces_of_v = defaultdict(set)
    for i, p in piece_of_edge.items():
        if p is not None:
            pieces_of_v[i].add(p)
            pieces_of_v[(i + 1) % L].add(p)
    for fv in faces:
        on = [pos[u] for u in fv if u in pos]
        fset = set(fv)
        # face edges on C'
        cedges = set()
        for a, b in zip(fv, fv[1:] + fv[:1]):
            if a in pos and b in pos:
                i, j = pos[a], pos[b]
                if (i + 1) % L == j:
                    cedges.add(i)
                elif (j + 1) % L == i:
                    cedges.add(j)
        Vs = [i for i in on if lab[i] == 'V']
        if len(Vs) >= 2:
            if not (len(Vs) == 2 and any(piece_of_edge[i] is None and {i, (i + 1) % L} == set(Vs) for i in cedges)):
                return False
        bypiece = defaultdict(set)
        for i in on:
            for p in pieces_of_v[i]:
                bypiece[p].add(i)
        for p, s in bypiece.items():
            if len(s) >= 3:
                return False
            if len(s) == 2:
                i, j = sorted(s)
                e = i if (i + 1) % L == j else (j if (j + 1) % L == i else None)
                if e is None or e not in cedges or piece_of_edge[e] != p:
                    return False
    return True


def double(H, link, lab):
    L = len(link)
    pos = {u: i for i, u in enumerate(link)}
    isV = {u: (u in pos and lab[pos[u]] == 'V') for u in H}
    M = nx.MultiGraph()

    def cp(u, side):
        return (u, 'M') if isV[u] else (u, side)

    for a, b in H.edges():
        mirror = isV[a] and isV[b] and a in pos and b in pos and (
            (pos[a] + 1) % L == pos[b] or (pos[b] + 1) % L == pos[a])
        if mirror:
            M.add_edge(cp(a, 'A'), cp(b, 'A'))
        else:
            M.add_edge(cp(a, 'A'), cp(b, 'A'))
            M.add_edge(cp(a, 'B'), cp(b, 'B'))
    for i, u in enumerate(link):
        if lab[i] == 'X':
            M.add_edge((u, 'A'), (u, 'B'))
    for u in H:
        M.add_node(cp(u, 'A'))
        M.add_node(cp(u, 'B'))
    return M


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    weights = defaultdict(Fraction)
    agree = disagree = 0
    for n in range(4, nmax + 1):
        for line in open(f'planar2c_{n}.g6', 'rb'):
            G = nx.from_graph6_bytes(line.strip())
            if nx.node_connectivity(G) < 3:
                continue
            darts, sig, alp = embed(G)
            faces = cells(darts, sig, alp)
            A = len(automorphisms(sig, alp))
            for w in G:
                fw = [f for f in faces if any(darts[d][0] == w for d in f)]
                if any(len(f) != 3 for f in fw):
                    continue
                link = link_cycle(darts, sig, w)
                H = G.copy()
                H.remove_node(w)
                inner = inner_faces_of_H(darts, faces, w)
                for lab in itertools.product('VXP', repeat=len(link)):
                    ok = conditions(lab, link, inner)
                    M = double(H, link, lab)
                    simple = not any(M.number_of_edges(a, b) > 1 for a, b in M.edges())
                    Ms = nx.Graph(M)
                    good = simple and Ms.number_of_nodes() >= 4 and nx.node_connectivity(Ms) >= 3
                    if good == ok:
                        agree += 1
                    else:
                        disagree += 1
                        if disagree <= 10:
                            print('DISAGREE', line.strip(), w, link, ''.join(lab), 'cond', ok, '3conn', good)
                    if good:
                        v = Ms.number_of_nodes()
                        f = Ms.number_of_edges() - v + 2
                        nv = sum(1 for c in lab if c == 'V')
                        nm = sum(1 for i in range(len(lab)) if lab[i] == 'V' and lab[(i + 1) % len(lab)] == 'V')
                        nc = sum(1 for c in lab if c == 'X')
                        weights[(v, f, nv, nm, nc)] += Fraction(1, 2 * A)
    print('agree', agree, 'disagree', disagree)
    stats = {}
    for fn in glob.glob('stats_*.json'):
        for key, d in json.load(open(fn)).items():
            stats.setdefault(key, {}).update(d)
    bad = ok = 0
    vmax = nmax - 1
    keys = set()
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v > vmax:
            continue
        for t, wgt in d.items():
            t = json.loads(t)
            if t[0] == 'ref':
                keys.add((v, f, t[1], t[2], t[3]))
    for (v, f, nv, nm, nc), wgt in weights.items():
        if v <= vmax:
            keys.add((v, f, nv, nm, nc))
    for k in sorted(keys):
        v, f, nv, nm, nc = k
        want = Fraction(stats.get(f'{v},{f}', {}).get(json.dumps(['ref', nv, nm, nc]), '0'))
        got = weights.get(k, 0)
        if want == got:
            ok += 1
        else:
            bad += 1
            print('weight mismatch', k, 'half-count', got, 'brute', want)
    print('weights ok', ok, 'bad', bad)


if __name__ == '__main__':
    main()
