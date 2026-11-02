# A076906 family: meanders crossing two perpendicular roads

A river (a simple curve from infinity to infinity) crosses the two axes of the plane n times, never through
the origin. Curves count up to isotopy that keeps the axes in place. Six OEIS entries (Sloane and Wild 2002–03,
partly extended by Robert Price 2012 to n = 20):

| entry | counts | before |
|---|---|---|
| A076906 | directed, from (−,−) to (+,+) for even n, to (+,−) for odd n | n ≤ 20 |
| A076907 | directed, from (−,−) to any quadrant | n ≤ 20 |
| A077551 | as A076907, first crossing on the x axis | n ≤ 20 |
| A076875 | undirected, one end in (−,−) | n ≤ 7 |
| A077550 | undirected, all quadrants, up to the 8 symmetries of the axes | n ≤ 7 |
| A077556 | undirected, one end in (−,−), first crossing (from a (−,−) end) on the x axis | n ≤ 7 |

**Result:** all six to **n = 37** (A077556 to n = 35) by a transfer matrix (minutes per run), plus exact formulas for A077551, A076875,
A077550 and A077556 in terms of the directed counts. Four old terms are wrong:

* A076875(6) = 176 should be **146** = A076906(6) + (A076907(6) − A076906(6))/2 = 82 + 64.
  The listed value contradicts the (corrected) terms of A076906/A076907.
* A077550(5) = 21 and A077550(7) = 131 should be **16** and **110** (Burnside below; also direct
  orbit enumeration, orbits.py).
* A077556(7) = 239 should be **218** = A077551(7) = A076906(7). 239 is probably the pre-2003 wrong value
  of A076906(7) carried over ("a(7) corrected Aug 25 2003").

Pictures (every river for n ≤ 8, a viewer with each river's itinerary, the table of counts):
https://guillochon.github.io/rivers-at-the-crossroads/

## 1. Combinatorial model (brute.py)

Cut the plane along the four rays W, S, E, N. Inside each quadrant the river pieces form a non-crossing
perfect matching of the quadrant's boundary points read in linear order: crossings on the first ray going
outward, the river ends at infinity (rim) in this quadrant, crossings on the second ray coming back in.
The corner at O carries no point, so cyclic order is linear. An isotopy class is exactly a choice of crossing
counts per ray and one matching per quadrant whose union is one path from start to end. Two ends in the same
quadrant come in two orders along the rim (so A076907(0) = 2). brute.py enumerates this directly; it
reproduces A076906 and A076907 for n ≤ 8.

## 2. Transfer matrix (tm.py prototype, rust/ engine)

The processed part of the plane lies left of a cut curve. The cut is a sequence (top to bottom) of markers,
where it crosses a ray, and items, the river edges crossing it. Two items joined inside the processed part
are a non-crossing pair (written as parentheses); START / END items lead to the river's ends. Stages, each
processing the points of one ray in a fixed direction:

1. cut x = c < 0, marker W moves right from −∞ to O: `[Q2 items] W [Q3 items]`
2. marker N enters at the top and moves down from +∞ to O: `[Q1] N [Q2] W [Q3]`
3. when N and W are adjacent at O (nothing between), push the cut across O: `N W` → `E S`;
   marker S moves down from O to −∞: `[Q1] E [Q4] S [Q3]`; at the end nothing may remain below S
4. cut x = c > 0, marker E moves right to +∞: `[Q1] E [Q4]`

At each crossing point the river has one half-edge on each side of the active marker, and each either closes
the item next to the marker on its side or opens a new one. This is Jensen's meander step. Closed loops are
forbidden. Pruning: Q1 items close at N or E points, Q2 at W/N, Q3 at W/S, Q4 at S/E, one item per side
per point, so e.g. in stage 2 the remaining points must be at least |Q3| + max(|Q2|, |Q1|).

State count grows about 2.5× per two crossings. n = 36 peaks at 40.5M states (4 min); n = 37 at about
60M. Memory (u128 key + u128 count) limits this 32 GB machine to roughly n = 38.

Runs: `run_all.sh 36` (ends in Q1 and Q3 for even n, in Q4 for odd n). Odd n need only Q4, since the
counts for Q4 and Q2 agree by reflection in y = x.

## 3. Formulas

Reflection in y = x fixes Q3 and Q1, swaps Q4 and Q2, and swaps the axes:

* A077551(n) = A076907(n)/2 for even n, = A076906(n) for odd n (all 21 known terms agree).
* A076875(n) = A076906(n) + (A076907(n) − A076906(n))/2 for even n, = 2·A076906(n) for odd n
  (undirected Q3→Q3 curves are counted twice in A076907).
* A077556(n) = A076906(n) for odd n. For even n, A077556(n) = A076906(n)/2 + (C3(n) − Cxx(n))/2, where C3
  counts directed Q3→Q3 curves and Cxx those whose first crossings from both ends are on the x axis
  (engine option `first`).

**A077550, by Burnside.** Rotations fix no curve: a fixed curve's midpoint would have to be O. An axis reflection
fixes a curve only for odd n: the midpoint lies on the axis and each half crosses only the other road. A
diagonal reflection fixes a curve only for even n. The halves are counted by the engine with rays switched
off (`noW noS noE`, `noS noE`), and they turn out to be known sequences:

* half curves for an axis reflection, k crossings: **A005316(k)** (open meanders);
* half curves for a diagonal reflection: **A076876(k)** (meanders across two parallel roads).

Both identities were checked against the OEIS b-files for k ≤ 20. They are also plain to see: an
upper half-plane cut by the N ray, entered from the x axis, is an open meander. The half-plane above y = x
holding the rays W and N is the strip picture of two parallel roads. Hence, for n ≥ 1,

    A077550(n) = (A076907(n) + 2·A005316((n−1)/2)) / 4    (n odd)
    A077550(n) = (A076907(n) + 2·A076876(n/2)) / 4        (n even)

## 4. Checks

* brute.py = engine for n ≤ 8; orbits.py (explicit curves, symmetry orbits) gives A076875, A077550, A077556
  for n ≤ 8 and agrees with the formulas.
* All 21 known terms of A076906, A076907, A077551 (Price) reproduced.
* Reflection in the x axis maps Q3→Q1 / Q3→Q4 / Q3→Q3 to Q2→Q4 / Q2→Q1 / Q2→Q2. The engine starting in Q2
  takes different routes through the stages (a START item above W, ends in other regions). Full runs agree
  term by term: results/q2*.txt.
* Half-curve counts match A005316 and A076876 (k ≤ 20).

## Files

brute.py, orbits.py (reference enumeration), tm.py (Python prototype), rust/ (engine `meander2`),
run_all.sh, combine.py (→ results/counts.json, bfiles/). The gallery page's source is in the private project;
its published copy is the rivers-at-the-crossroads repo.
