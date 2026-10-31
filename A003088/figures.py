"""Figures for the A003088 write-up (matplotlib).

    python figures.py     # writes a003088_n3.png (all digraphs on 3 nodes) and a003088_growth.png (from b003088.txt)
"""
import itertools
from math import factorial
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

HERE = Path(__file__).resolve().parent
N = 3
PAIRS = [(i, j) for i in range(N) for j in range(N) if i != j]
POS = [(0, 0.95), (-0.85, -0.5), (0.85, -0.5)]
COLS = ["#2a78d6", "#eb6834", "#1baf7a"]


def reach(arcs):
    R = [[i == j for j in range(N)] for i in range(N)]
    for i, j in arcs:
        R[i][j] = True
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if R[i][k] and R[k][j]:
                    R[i][j] = True
    return R


def canon(arcs):
    return min(tuple(sorted((p[i], p[j]) for i, j in arcs)) for p in itertools.permutations(range(N)))


def classify(arcs):
    """(unilateral?, strong components in chain order)"""
    R = reach(arcs)
    uni = all(R[i][j] or R[j][i] for i in range(N) for j in range(N))
    comps = []
    for v in range(N):
        for c in comps:
            if R[v][c[0]] and R[c[0]][v]:
                c.append(v)
                break
        else:
            comps.append([v])
    comps.sort(key=lambda c: -sum(R[c[0]]))  # earlier in the chain reaches more nodes
    return uni, comps


def draw(ax, arcs, ok, comps):
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.05, 1.45); ax.set_aspect("equal"); ax.axis("off")
    cid = {v: k for k, c in enumerate(comps) for v in c}
    s = set(arcs)
    for i, j in arcs:
        rad = 0.25 if (j, i) in s else 0.0
        ax.add_patch(FancyArrowPatch(POS[i], POS[j], connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                                     mutation_scale=13, lw=1.6, color="#3a3a38" if ok else "#b4b2a9",
                                     shrinkA=13, shrinkB=13))
    for v in range(N):
        ax.add_patch(Circle(POS[v], 0.2, color=COLS[cid[v]] if ok else "#d3d1c7", zorder=3))
        if ok and len(comps) > 1:
            ax.text(*POS[v], str(cid[v] + 1), ha="center", va="center", color="white", fontsize=9,
                    fontweight="bold", zorder=4)


def figure_n3():
    graphs = {}
    for mask in range(1 << len(PAIRS)):
        c = canon([PAIRS[k] for k in range(len(PAIRS)) if mask >> k & 1])
        graphs.setdefault(c, list(c))
    data = [(g, *classify(g)) for g in graphs.values()]
    uni = sorted([d for d in data if d[1]], key=lambda d: (len(d[0]), -len(d[2])))
    non = sorted([d for d in data if not d[1]], key=lambda d: len(d[0]))
    assert len(uni) == 11 and len(non) == 5  # a(3) = 11 of the 16 digraphs on 3 nodes

    fig = plt.figure(figsize=(7.2, 9.6))
    fig.suptitle("What A003088 counts: unilateral digraphs, n = 3", fontsize=13, y=0.985)
    fig.text(0.5, 0.945, "For every two nodes, at least one can reach the other along arrows.\n"
             "Colors = strong components; numbers = their order in the chain 1 → 2 → 3.",
             ha="center", va="top", fontsize=9.5, color="#52514e")
    fig.text(0.04, 0.87, "a(3) = 11 unilateral", fontsize=11, fontweight="bold", color="#1b2330")
    for k, (g, ok, comps) in enumerate(uni):
        r, c = divmod(k, 4)
        ax = fig.add_axes([0.03 + c * 0.24, 0.68 - r * 0.18, 0.22, 0.17])
        draw(ax, g, ok, comps)
        ax.set_title(f"{len(g)} arcs" + ("  (strong)" if len(comps) == 1 else ""), fontsize=8.5, color="#52514e", pad=0)
    fig.text(0.04, 0.275, "Not counted: 5 digraphs where two nodes cannot reach each other either way",
             fontsize=10.5, fontweight="bold", color="#898781")
    for k, (g, ok, comps) in enumerate(non):
        ax = fig.add_axes([0.01 + k * 0.196, 0.08, 0.18, 0.17])
        draw(ax, g, ok, comps)
        ax.set_title(f"{len(g)} arcs", fontsize=8.5, color="#898781", pad=0)
    fig.text(0.5, 0.025, "Every automorphism keeps each strong component in place, since the chain order is unique.\n"
             "That is what lets the formula count these up to n = 50 without listing them.",
             ha="center", fontsize=9, color="#52514e")
    fig.savefig(HERE / "a003088_n3.png", dpi=170)


def figure_growth():
    a = [int(line.split()[1]) for line in open(HERE / "b003088.txt")]
    dig = [len(str(x)) for x in a]
    n = list(range(len(a)))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7.5), gridspec_kw={"height_ratios": [3, 2]})
    for ax in (ax1, ax2):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color="#e1e0d9"); ax.set_axisbelow(True)
    ax1.plot(n[7:], dig[7:], color="#2a78d6", lw=2, marker="o", ms=3, label=f"new, by formula (n = 8–{n[-1]})")
    ax1.plot(n[:8], dig[:8], "o", color="#898781", ms=6, label="known before (n ≤ 7)")
    ax1.set_xlabel("n (nodes)"); ax1.set_ylabel("digits of a(n)")
    ax1.set_title("A003088: unilateral digraphs on n unlabeled nodes", fontsize=12)
    ax1.annotate(f"a({n[-1]}): {dig[-1]} digits", (n[-1], dig[-1]), xytext=(28, 0.89 * dig[-1]),
                 arrowprops=dict(arrowstyle="->", color="#52514e"), fontsize=10)
    ax1.legend(frameon=False, loc="upper left")
    m = range(3, 15)
    r = [a[k] * factorial(k) / 2 ** (k * (k - 1)) for k in m]
    ax2.axhline(1, color="#c3c2b7", lw=1)
    ax2.plot(list(m), r, color="#2a78d6", lw=2, marker="o", ms=4)
    ax2.set_xlabel("n"); ax2.set_ylabel("a(n) · n! / 2^(n(n−1))")
    ax2.set_title("Tends to 1: almost every digraph is strong", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "a003088_growth.png", dpi=160)


if __name__ == "__main__":
    figure_n3()
    figure_growth()
