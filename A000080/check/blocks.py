"""Check the block structure behind the formula on explicit lists of minimal triangle graphs.

    ./a000080 -t n > t{n}.txt          # Martin Fuller's program (check/a000080_fuller.c), n = 3..15
    python check/blocks.py DIR NMAX    # DIR holds t3.txt .. t{NMAX}.txt

Splits every graph into the blocks of its point-triangle incidence graph (triangles in the same 2-connected
piece go together), keeps one copy of each block up to isomorphism, and prints the number of distinct blocks by
number of points.  With NMAX = 15 this gives 1, 1, 1, 2, 2, 4, 6, 12, 20, 42, 84, 187, 424 for 3..15 points,
which is what species.py / rust predict for B*(x, x): one-edge blocks plus 2-connected loopless multigraphs,
counted by vertices + edges.
"""
import re
import sys
from collections import Counter

import networkx as nx
from networkx.algorithms import isomorphism as iso


def blocks(T):
    G = nx.Graph()
    for i, t in enumerate(T):
        for v in t:
            G.add_edge(('t', i), ('v', v))
    uf = list(range(len(T)))

    def f(x):
        while uf[x] != x:
            x = uf[x]
        return x

    for comp in nx.biconnected_components(G):
        ts = [x[1] for x in comp if x[0] == 't']
        for a in ts[1:]:
            uf[f(a)] = f(ts[0])
    groups = {}
    for i in range(len(T)):
        groups.setdefault(f(i), []).append(T[i])
    return list(groups.values())


def incidence(B):
    G = nx.Graph()
    for i, t in enumerate(B):
        G.add_node(('t', i), k=1)
        for v in t:
            G.add_node(('v', v), k=0)
            G.add_edge(('t', i), ('v', v))
    return G


def main(d, nmax):
    nm = iso.categorical_node_match('k', 0)
    reps = {}
    for n in range(3, nmax + 1):
        for line in open(f'{d}/t{n}.txt').read().split():
            T = [tuple(map(int, x)) for x in re.findall(r'\((\d+),(\d+),(\d+)\)', line)]
            for B in blocks(T):
                G = incidence(B)
                k = len({v for t in B for v in t})
                h = nx.weisfeiler_lehman_graph_hash(G, node_attr='k')
                lst = reps.setdefault((k, h), [])
                if not any(nx.is_isomorphic(G, H, node_match=nm) for H in lst):
                    lst.append(G)
    c = Counter()
    for (k, _), lst in reps.items():
        c[k] += len(lst)
    print(sorted(c.items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
