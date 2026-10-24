"""Plot a(n) (GG_n pieces on the plateau above sec(pi/n)) against the period-4 cubic fits, log-log."""
import matplotlib.pyplot as plt
import numpy as np

FIT = {0: (1, 19, -32, 40), 1: (1, 73, -145, 127), 2: (1, -1, 0, 4), 3: (1, 77, -93, 79)}
COLORS = {0: "#2a78d6", 1: "#eb6834", 2: "#1baf7a", 3: "#eda100"}   # dataviz categorical slots 1-4
MARKERS = {0: "o", 1: "s", 2: "^", 3: "D"}
INK, MUTED, GRID, SURFACE = "#1f1f1e", "#6b6a64", "#e6e5e0", "#fcfcfb"
ARTIFACT = set()   # none since the exact solver (ExactPuzzle)


a = {int(n): int(v) for n, v in (l.split() for l in open("results/a_n.txt"))}

fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, height_ratios=(3, 1.4),
                              facecolor=SURFACE, layout="constrained")
for axis in (ax, ax2):
    axis.set_facecolor(SURFACE)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(MUTED)
    axis.tick_params(colors=MUTED, labelcolor=INK)
    axis.grid(True, which="major", color=GRID, linewidth=0.8)
    axis.set_axisbelow(True)

xs = np.linspace(4, 210, 600)
ax.plot(xs, xs**3 / 8, color=MUTED, lw=1.2, ls=(0, (4, 3)), label="n³/8", zorder=1)
for c in range(4):
    ns = sorted(n for n in a if n % 4 == c and n >= 4)
    ax.plot(xs, np.polyval(FIT[c], xs) / 8, color=COLORS[c], lw=2, alpha=0.55, zorder=2)
    good = [n for n in ns if n not in ARTIFACT]
    ax.scatter(good, [a[n] for n in good], s=34, marker=MARKERS[c], color=COLORS[c], edgecolor=SURFACE,
               linewidth=1.2, zorder=3, label=f"n ≡ {c} (mod 4)")
    bad = [n for n in ns if n in ARTIFACT]
    if bad:
        ax.scatter(bad, [a[n] for n in bad], s=60, marker=MARKERS[c], facecolor="none", edgecolor=INK,
                   linewidth=1.2, zorder=4)
    # residual panel: data relative to n^3/8
    ax2.plot(xs, np.polyval(FIT[c], xs) / xs**3, color=COLORS[c], lw=2, alpha=0.55)
    ax2.scatter(good, [8 * a[n] / n**3 for n in good], s=22, marker=MARKERS[c], color=COLORS[c],
                edgecolor=SURFACE, linewidth=1, zorder=3)

ax.scatter([3], [a[3]], s=60, marker="*", color=INK, zorder=4)
ax.annotate("n = 3 (r = 2, degenerate)", (3, a[3]), xytext=(8, 6), textcoords="offset points", color=MUTED,
            fontsize=9)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylabel("pieces a(n)", color=INK)
ax.set_title("Two-disk circle puzzle GG_n: pieces just past r = sec(π/n)", color=INK, loc="left", fontsize=12, pad=20)
ax.text(0, 1.012, "points = computed (n = 3..200), lines = period-4 cubic fits", transform=ax.transAxes,
        color=MUTED, fontsize=9, va="bottom")
ax.legend(frameon=False, loc="upper left", fontsize=9, labelcolor=INK)
ax2.axhline(1, color=MUTED, lw=1, ls=(0, (4, 3)))
ax2.set_ylabel("8·a(n) / n³", color=INK)
ax2.set_xlabel("n (rotation order)", color=INK)
ax2.set_xlim(2.8, 240)
ax2.set_xticks([3, 5, 10, 20, 50, 100, 200])
ax2.set_xticklabels(["3", "5", "10", "20", "50", "100", "200"])
fig.savefig("results/an_loglog.png", dpi=160, facecolor=SURFACE)
print("saved", len(a), "terms")
