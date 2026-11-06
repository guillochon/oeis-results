"""Brute-force two-sort cycle index of 3-connected planar graphs (pilot input for the tower).

For every 3-connected planar graph K (read from graph6 files of 2-connected planar graphs with
min degree 3, made by ../../A000944/brute/gen_poly.sh) and every automorphism g of K (by Whitney,
the map automorphisms of its embedding, reflections included) record the monomial
    (vertex cycle type of g,  edge orbits of g as (length k, flipped?))
where an edge orbit of length k is "flipped" when g^k reverses each of its edges.
Output: JSON list of [n, vertex_partition, edge_orbits, weight] with weight = sum of 1/|Aut K|,
so summing weights over all monomials of a given n gives A000944(n).

Usage: python zk.py OUT.json file1.g6 [file2.g6 ...]
"""
import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from multiprocessing import Pool

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
                out.append(m)
    return out


def monomial(m, darts, alp, nverts):
    vert_of = [u for (u, w) in darts]
    # vertex permutation
    vp = {}
    for i in range(len(m)):
        vp[vert_of[i]] = vert_of[m[i]]
    seen = set()
    vparts = []
    for v in range(nverts):
        if v in seen:
            continue
        k = 0
        x = v
        while x not in seen:
            seen.add(x)
            x = vp[x]
            k += 1
        vparts.append(k)
    # edge orbits
    seen = set()
    eorb = []
    for d in range(len(m)):
        if d in seen:
            continue
        k = 0
        x = d
        while True:
            seen.add(x)
            seen.add(alp[x])
            x = m[x]
            k += 1
            if x == d or x == alp[d]:
                break
        eorb.append((k, int(x == alp[d])))
    return tuple(sorted(vparts, reverse=True)), tuple(sorted(eorb))


def work(lines):
    acc = defaultdict(Fraction)
    for line in lines:
        G = nx.from_graph6_bytes(line)
        if nx.node_connectivity(G) < 3:
            continue
        G = nx.convert_node_labels_to_integers(G)
        darts, sig, alp = embed(G)
        auts = automorphisms(sig, alp)
        A = len(auts)
        n = G.number_of_nodes()
        for m in auts:
            vp, eo = monomial(m, darts, alp, n)
            acc[(n, vp, eo)] += Fraction(1, A)
    return acc


def main():
    try:
        import psutil
        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass
    out = sys.argv[1]
    lines = []
    for fn in sys.argv[2:]:
        with open(fn, 'rb') as fh:
            lines += [l.strip() for l in fh if l.strip()]
    chunks = [lines[i:i + 2000] for i in range(0, len(lines), 2000)]
    tot = defaultdict(Fraction)
    with Pool(6) as pool:
        for acc in pool.imap_unordered(work, chunks):
            for k, w in acc.items():
                tot[k] += w
    rows = [[n, list(vp), [list(e) for e in eo], str(w)] for (n, vp, eo), w in sorted(tot.items())]
    with open(out, 'w') as fh:
        json.dump(rows, fh)
    byn = defaultdict(Fraction)
    for (n, vp, eo), w in tot.items():
        byn[n] += w
    print({n: str(w) for n, w in sorted(byn.items())})


if __name__ == '__main__':
    main()
