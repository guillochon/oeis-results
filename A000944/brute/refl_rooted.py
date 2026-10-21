"""Brute-force anchor-rooted counts of valid labelled halves, by (v_M, f_M, #anchors) and by the
root node type, to debug formula/refl.py.  rooted = sum_W sum_{w, lab valid} 2A / |Aut W|.

Usage: python refl_rooted.py NMAX_CONE > refl_rooted_N.json
"""
import itertools
import json
import sys
from collections import defaultdict
from fractions import Fraction

import networkx as nx

from refl_check import conditions, double, inner_faces_of_H, link_cycle
from symstats import automorphisms, cells, embed


def main():
    nmax = int(sys.argv[1])
    out = defaultdict(Fraction)
    for n in range(4, nmax + 1):
        for line in open(f'planar2c_{n}.g6', 'rb'):
            G = nx.from_graph6_bytes(line.strip())
            if nx.node_connectivity(G) < 3:
                continue
            darts, sig, alp = embed(G)
            faces = cells(darts, sig, alp)
            Aw = len(automorphisms(sig, alp))
            for w in G:
                fw = [f for f in faces if any(darts[d][0] == w for d in f)]
                if any(len(f) != 3 for f in fw):
                    continue
                link = link_cycle(darts, sig, w)
                H = G.copy()
                H.remove_node(w)
                inner = inner_faces_of_H(darts, faces, w)
                for lab in itertools.product('VXP', repeat=len(link)):
                    if not conditions(lab, link, inner):
                        continue
                    M = nx.Graph(double(H, link, lab))
                    v = M.number_of_nodes()
                    f = M.number_of_edges() - v + 2
                    A = sum(1 for c in lab if c in 'VX')
                    out[f'{v},{f},{A}'] += Fraction(2 * A, Aw)
    json.dump({k: str(v) for k, v in sorted(out.items())}, sys.stdout, indent=0)


if __name__ == '__main__':
    main()
