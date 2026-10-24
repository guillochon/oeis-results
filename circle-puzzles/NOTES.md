# Two-disk circle puzzles GG_n(r): piece counts

Setting of Hearn, Kretschmer, Rokicki, Streeter, Vergo, "Two-Disk Compound Symmetry Groups"
([arXiv:2302.12950](https://arxiv.org/abs/2302.12950)): disks of radius r centred at (-1, 0) and (1, 0), each turning
by 2*pi/n. The cut set is the closure of both boundary circles under the two rotations. Pieces are the bounded
faces of that arrangement, inside the union of the disks.

`twodisk.py` stores every cut as an arc of a radius-r circle (a centre plus a union of angle intervals) and closes
the set with a worklist. It counts pieces with Euler's formula (F = E - V + C). `Puzzle.check()` independently
traces every face and checks that exactly one face is unbounded and that the bounded faces' areas add up to the
area of the disk union. The two agree, to within 1e-15 in area, on every case tried.

- `staircase.py n lo hi m`: every change in the piece count for r in [lo, hi]. Changes closer together than 1e-5
  are merged, because slivers thinner than the point tolerance (1e-7) appear just past degenerate radii.
- `bisect_tr.py n lo hi`: finds one transition radius.
- `plateau2.py n0 n1`: a(n) below.

## Results

- r < 1: 2 pieces. 1 < r < sec(pi/n): 2n+1 pieces (the lens and its n-1 images in each disk). This matches the
  paper's "9 pieces added" for GG_5.
- **a(n) = pieces for r just above sec(pi/n)**, where neighbouring lens images first overlap. Each value is
  constant at r = sec(pi/n)*(1+1e-4) and *(1+1e-3), and for n <= 16 also at *(1+3e-3). For n = 3..48:

  69, 35, 169, 23, 443, 189, 683, 113, 1213, 515, 1597, 319, 2423, 1061, 2959, 689, 4121, 1875, 4817, 1271, 6355,
  3005, 7219, 2113, 9173, 4499, 10213, 3263, 12623, 6405, 13847, 4769, 16753, 8771, 18169, 6679, 21611, 11645,
  23227, 9041, 27245, 15075, 29069, 11903, 33703, 19109

  For n >= 4 it is a cubic quasi-polynomial with period 4. Each cubic is fitted on 4 terms and confirmed on the
  other 7-8 in its class:

  | n mod 4 | a(n) |
  |---|---|
  | 0 | (n^3 + 19n^2 - 32n + 40)/8 |
  | 1 | (n^3 + 73n^2 - 145n + 127)/8 |
  | 2 | (n^3 - n^2 + 4)/8 |
  | 3 | (n^3 + 77n^2 - 93n + 79)/8 |

  The plateau is wide (`plateau_width.py`): a(n) holds on (sec(pi/n), R_n] with R_n = 2, 2cos(pi/10), sqrt 3,
  2cos(3pi/14), sqrt 2, sec(2pi/9), sec(pi/5), sec(2pi/11), sec(pi/6) for n = 4..12 (sec(2pi/n) for n >= 8 so far),
  about 3-4 times sec(pi/n) - 1. 300-400 random radii inside each plateau (n = 5, 7, 8, 10, 12; `plateau_dense.py`)
  all gave a(n).

  n = 3 (69) is off the n = 3 (mod 4) cubic. There sec(pi/3) = 2, the radius at which each disk reaches the other's
  centre. The formulas are conjectural: they fit the data, but there is no proof.
- n = 5 staircase: 11 -> 169 at sqrt(5)-1 = sec(36 deg), 169 -> 269 at 2 sin 72 deg = 1.9021130326 (A188593),
  269 -> 287 at 1.9077001924, ... then growth toward the critical radius sqrt(3+phi) = 2.148961 from the paper
  (results/n5_crit.txt).
- My runs reached their circle cap before each paper critical radius (n=5: cap at about 2.075 < 2.149; n=7: 1.6025 <
  1.6236), so there is no conflict with the paper. The cap is a resource limit, not a critical radius.

## Larger n (pieces_fast)

`pieces_fast` replaces the O(N^2) Python loop over circle pairs: arcs are cut into sub-arcs of span <= 0.05, a k-d
tree on their bounding discs gives candidate pairs, and intersections, vertex clustering (k-d tree + connected
components) and edge counting are numpy. It is identical to `pieces()` on every comparison case, and n = 100 takes
about 1 s. `plateau3.py` samples r = sec + f*(sec - 1) for f = 1.3, 1.0, 1.6 (half the cores at most; use 2
workers, memory grows to ~1 GB per worker near n = 200) and face-traces the first radius.

- n = 3..200 (float solver): the three radii agree and match the cubic fits except n = 189, 195 (results/plateau3_*.txt,
  plot results/an_loglog.png from `plot_an.py`).
- Double precision is the limit. Near the bottom of the plateau (f = 0.2) odd n >= 85 have pieces of area ~5e-17,
  and face tracing fails there; at f = 1 odd n >= 87 still have 1e-12..1e-15 slivers, so the area check is moved to
  f = 1.3. Beyond n ~ 150 the area check fails broadly while the counts still agree across radii.
- Tolerances were tightened to EPS_PT = 1e-10, EPS_ANG = 1e-13 (no change for n <= 60). At n = 189 that is too
  tight: accumulated rotation error leaves near-duplicate arcs (208868 arcs vs 185180), giving 1161564; with
  (1e-9, 1e-12) the count is 1166453, exactly the fit.

## Exact closure (ExactPuzzle) and the final table

The drift came from computing each centre and arc end by rotating an already-rotated one. `ExactPuzzle` removes it:

- Centres lie in Z[zeta_n]; they are stored as integer vectors modulo the cyclotomic polynomial and identified by
  key, with no tolerance. Coefficients stay tiny (max 16 at n = 195).
- Arc ends are labels: ('c', p, s) is the crossing with circle p (sign s picks one of the two), ('f', k) is the
  angle 2*pi*k/n at the 0/2pi seam. A rotation maps labels to labels, and angles are always recomputed from exact
  centres, so nothing accumulates.

It reproduces every earlier reference count. n = 189 and 195 now give 1166453 and 1290593 (the fits), unchanged
for point tolerances 1e-9, 1e-10 and 1e-11, with 141751 arcs at n = 189 (no near-duplicates). The full run
(results/exact_3_200.txt, 3 workers) gives the same count at f = 1.3, 1.0 and 1.6 for every n = 5..200, every one
equal to the cubic fit, and the face-tracing area check now passes for every n (it failed for most odd n >= 87
before). `an_table.py` writes results/a_n.txt (n = 3..200, b-file style; n = 3, 4 from the sec*(1+1e-4) runs, as
their plateaus are narrower than f = 1.3) and fills the list in the viz page.

The page `viz/circle_puzzles/index.html` runs a JavaScript port of `ExactPuzzle` in the browser (identical circle
and arc counts; n = 200 in about 1.5 s) and draws the cut set for any n <= 200 with deep zoom.

## OEIS check (local dump of 2026-10-16)

- a(n) itself, each mod-4 class, and the transforms (a-1)/2, a-2 and a-2n-1 (DATA and b-file terms): no match.
- Critical radii sqrt(3+phi), sqrt(4-phi), sqrt(5(2-sqrt 2)) and sqrt(2(20-11 sqrt 3)) as decimal expansions: not
  in OEIS.
- No OEIS entry mentions circle puzzles, jumbling or Hearn fractals.
