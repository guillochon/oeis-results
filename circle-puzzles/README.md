# Pieces of the two-disk circle puzzle GG_n just past r = sec(π/n)

A new sequence, n = 3..200:

    69, 35, 169, 23, 443, 189, 683, 113, 1213, 515, 1597, 319, 2423, 1061, 2959, 689, 4121, 1875, ...

(all 198 terms in [`b_circle.txt`](b_circle.txt) and [`results/a_n.txt`](results/a_n.txt)). An interactive viewer
that recomputes every cut and piece in the browser, for any n ≤ 200, is at
**https://guillochon.github.io/circle-puzzle-cuts/**.

## 1. What is counted

Take two disks of radius r centred at (−1, 0) and (1, 0), and let each turn by 2π/n about its centre. This is the
setting of Hearn, Kretschmer, Rokicki, Streeter and Vergo, *Two-Disk Compound Symmetry Groups*
([arXiv:2302.12950](https://arxiv.org/abs/2302.12950)), who call the group of moves GG_n(r).

The **cut set** is the smallest set of curves that contains both disk boundaries and is carried into itself by
each turn: whenever a turn of one disk moves part of a cut, the moved copy is a cut too. Every cut is an arc of a
circle of radius r. The **pieces** are the regions into which the cut set divides the union of the two disks. For
n = 2, 3, 4, 6 the cut set is finite at every radius; for every other n it is finite below a critical radius and a
fractal above it (Hearn et al.).

For 1 < r < sec(π/n) the lens shared by the disks is copied n − 1 times around each disk without overlap, and there
are 2n + 1 pieces. At r = sec(π/n) neighbouring copies of the lens start to overlap, and the count jumps. It is then
**constant** on a whole interval of radii. **a(n) is the number of pieces for r just above sec(π/n).**

The plateau is wide. Measured endpoints (bisection to 10⁻⁹):

| n | a(n) | plateau (sec(π/n), R_n] | R_n | next count |
|---|---|---|---|---|
| 4 | 35 | (1.414214, 2] | 2 | 67 |
| 5 | 169 | (1.236068, 1.902113] | 2 cos(π/10) | 269 |
| 6 | 23 | (1.154701, 1.732051] | √3 | 53 |
| 7 | 443 | (1.109916, 1.563663] | 2 cos(3π/14) | 615 |
| 8 | 189 | (1.082392, 1.414214] | sec(2π/8) | 303 |
| 9 | 683 | (1.064178, 1.305407] | sec(2π/9) | 1230 |
| 10 | 113 | (1.051462, 1.236068] | sec(2π/10) | 626 |
| 11 | 1213 | (1.042217, 1.188702] | sec(2π/11) | 2570 |
| 12 | 515 | (1.035276, 1.154701] | sec(2π/12) | 559 |
| 13 | 1597 | (1.029928, 1.129362] | sec(2π/13) | 3418 |
| 14 | 319 | (1.025717, 1.109916] | sec(2π/14) | 1439 |
| 15 | 2423 | (1.022341, 1.094636] | sec(2π/15) | 3829 |
| 16 | 1061 | (1.019591, 1.082392] | sec(2π/16) | 1668 |

So for 8 ≤ n ≤ 16 the plateau is exactly sec(π/n) < r ≤ sec(2π/n). For n = 3 the plateau starts at sec(π/3) = 2,
where each disk's boundary passes through the other's centre; a(3) = 69 holds just above 2 and the count changes
soon after.

## 2. Formulas (conjectured)

For n ≥ 4 every term fits a cubic quasi-polynomial of period 4:

| n mod 4 | a(n) |
|---|---|
| 0 | (n³ + 19n² − 32n + 40)/8 |
| 1 | (n³ + 73n² − 145n + 127)/8 |
| 2 | (n³ − n² + 4)/8 |
| 3 | (n³ + 77n² − 93n + 79)/8 |

Each cubic was fitted to the first 4 terms of its class and then checked against the other 45 to 46; all 195 terms
n = 4..200 agree. Equivalently, a(n) = a(n−1) + 3a(n−4) − 3a(n−5) − 3a(n−8) + 3a(n−9) + a(n−12) − a(n−13) for
n ≥ 17, and the generating function (offset 3, including the exceptional a(3)) is

    (69 − 34x + 134x² − 146x³ + 213x⁴ − 152x⁵ + 92x⁶ − 132x⁷ + 47x⁸ − 38x⁹ + 2x¹⁰ − 6x¹¹ − 5x¹² + 4x¹³)
      / ((1 − x)(1 − x⁴)³),

which reproduces all 198 terms. None of this is proved. a(n) ~ n³/8.

## 3. How the terms were computed

`twodisk.py`, class `ExactPuzzle`.

- **Exact closure.** Every cut is an arc of a radius-r circle whose centre is an image of ±1 under the turns, so the
  centres lie in Z[ζ_n]. They are stored as integer vectors modulo the n-th cyclotomic polynomial and identified by
  value, not by a floating-point tolerance (coefficients stay tiny: at most 16 at n = 195). An arc's ends are stored
  symbolically, as "where this circle crosses that circle" (or the angle 2πk/n at the 0/2π seam), so a turn maps
  ends to ends and every angle is recomputed from exact centres. Nothing drifts however long the orbit. A worklist
  applies each disk's turn to every new piece of cut inside that disk until nothing new appears.
- **Counting pieces.** The arcs are cut at all their crossings (candidate pairs from a k-d tree of short sub-arcs),
  coincident points are joined, and the pieces are the bounded faces of the resulting planar graph: F = E − V + C.
- **Radius.** Each term is computed at r = sec(π/n) + f·(sec(π/n) − 1) for f = 1.3, 1.0 and 1.6 (n ≥ 5), all
  inside the plateau, and the three counts must agree. n = 3, 4 have narrower plateaus and use r = sec(π/n)·(1 + 10⁻⁴)
  and (1 + 10⁻³) (also checked at 1 + 3·10⁻³).

## 4. Checks

- **Face tracing.** At the first radius every face of the arrangement is traced (half-edges sorted around each
  vertex by direction and curvature), and the run is accepted only if there is exactly one unbounded face, the
  number of bounded faces equals the Euler count, and the bounded faces' areas (Green's theorem on the arcs) sum to
  the area of the union of the two disks to within 10⁻⁹ (in practice 10⁻¹² or better). Every n = 5..200 passes
  (`results/exact_3_200.txt`, column "check-ok").
- **Three radii** per n give the same count (above).
- **Plateau scans.** For n = 4..16 the count was followed across the whole plateau by a grid plus bisection (table
  above), and 300 to 400 random radii inside the plateau for n = 5, 7, 8, 10, 12 all gave a(n).
- **An independent implementation.** The viewer is a JavaScript port written separately for the browser (crossings
  on a grid instead of a k-d tree, its own vertex hashing and face tracing). It reproduces the Python circle and arc
  counts exactly and counts 1,094,205 pieces at n = 200, equal to a(200).
- **Reference cases.** 2n + 1 pieces below sec(π/n), matching the paper's "9 pieces added" for GG_5; counts at
  other radii (e.g. 26,651 pieces for GG_7(1.59), 17,543 for GG_5(2.03)) agree between the two implementations.
- **Why exact arithmetic.** An earlier floating-point version (kept in `twodisk.py` as `Puzzle`) agreed for
  n ≤ 188 but went wrong at n = 189 and 195: accumulated rotation error left near-duplicate arcs. `ExactPuzzle`
  removed the discrepancy; both values now match the formulas.

## 5. Reproduce

Python 3 with numpy, scipy and sympy.

    python twodisk.py 7 1.2                 # quick check: (443, ...) pieces for n = 7 on the plateau (float solver)
    python plateau3.py 3 200 3              # all terms, 3 worker processes (~1 h; ~1 GB per worker near n = 200)
    python an_table.py                      # results/exact_3_200.txt -> results/a_n.txt
    python plot_an.py                       # results/an_loglog.png

From Python, `ExactPuzzle(n, r).close()` then `.pieces_fast()` gives one count and `.check()` the face-tracing
check.

## Files

- `twodisk.py`: the closure (`ExactPuzzle`; the float `Puzzle` it replaced), piece counting, face tracing.
- `plateau3.py`, `pool.py`: the run over n (half the cores, below-normal priority).
- `an_table.py`, `plot_an.py`: the table and the plot.
- `results/exact_3_200.txt`: the run log (three counts, the fit, the face check, time per n).
- `results/a_n.txt`, `b_circle.txt`: the terms.
- `NOTES.md`: the working log, including the plateau scans and the precision problems on the way.

Made by James Guillochon with Claude.
