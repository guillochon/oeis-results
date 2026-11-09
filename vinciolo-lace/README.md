# Symmetric lace designs without pinholes, and the minimum-region table

This folder computes two new OEIS entries. The first is a(1)–a(47) of a sequence:

> **a(n)** = number of n × n binary arrays that are invariant under the 8 symmetries of the square, in which every
> 1 is adjacent to another 1 horizontally, vertically or diagonally.

1, 2, 6, 7, 48, 55, 705, 861, 21454, 26696, 1311647, 1650868, 159700301, 203876423, … (`b_pinhole_d4.txt` has
n = 1..47). Without the symmetry, the same condition gives A202788.

The second is a table that contains it as column 2 (sections 6–8):

> **T(n,k)** = number of n × n binary arrays invariant under the 8 symmetries of the square in which every
> component of 1's (horizontal, vertical and diagonal adjacency) has at least k cells.

It is read by antidiagonals, T(1,1), T(1,2), T(2,1), …; `b_minregion_table.txt` has the antidiagonals n + k = 2..34
(561 terms).

It also contains the study of Federico de Vinciolo's lace patterns that led to the sequence. Pictures and
background are on the companion page, [Vinciolo's Lace](https://guillochon.github.io/vinciolo-lace/).

## 1. Where the rule comes from

Vinciolo's *Les Secondes Oeuvres, et Subtiles Inventions De Lingerie* (Paris, Jean Le Clerc, 1603) is a woodcut
pattern book. Most of its plates are charts for *lacis* (filet lace): a square knotted net in which chosen meshes
are darned solid. On the plates a printed square is a darned mesh and a blank square is an open one. Many plates,
headed "Carrez contenant 41 maille", hold two square panels each, drawn so that they look the same after any
rotation or reflection of the square.

The Metropolitan Museum of Art has the whole book online as public-domain photographs (accession 39.30). We traced
every plate with square geometric designs and kept the panels that are fully symmetric:

* **Tracing** (`tools/trace.py`). Dark squares are found as connected blobs; each is placed on an integer lattice by
  walking between neighboring blobs, rounding each step to the local cell pitch, which follows the curve of the page.
  Cells with no blob are open. Fourteen of the 41-mesh plates traced with no lattice conflicts and two with a few;
  the plates with a finer mesh (pages 17–25) and the checkered plates (42–47) are not used.
* **Panels** (`tools/panels.py`). A panel is a square window whose outer ring lies on Vinciolo's solid frame line
  (at least 90% printed) and which agrees with its 7 rotated and reflected copies on at least 90% of cells (83% for
  the upper panel of page 37). This finds 19 panels from 43 × 43 to 61 × 61 meshes on 13 pages. Each is
  symmetrized by a majority vote over its 8 copies, which removes tracing errors (and Vinciolo's rare slips).
  `panels.txt` lists all 19.
* **Rules** (`tools/rules.py`, `rare.py`, `props.py`, `loo.py`). With the frame counted as darned, the 19 panels
  contain every one of the 512 possible 3 × 3 blocks except one: an open mesh whose 8 neighbors are all darned.
  Random symmetric designs of the same density would contain about 120 such pinholes in the same area; the panels
  contain none. The next rarest block occurs 16 times, in one panel. Other candidate rules fail: the open meshes
  form one connected region in only 10 of the 19 panels, isolated darned dots occur in 7 (they are flower centers),
  and "every 3 × 4 block is one Vinciolo uses" is broken by 18 of the 19 panels when that panel is left out of the
  list. So the rule the panels never break is: **no open mesh is enclosed on all eight sides by darned meshes**.

Writing 1 for an open mesh and 0 for a darned one, with the frame darned (0 outside the array), the rule says
exactly that every 1 has a neighboring 1. Counting the symmetric designs that obey it gives a(n).

## 2. Reducing to one eighth of the array

Let h = ⌈n/2⌉ and fold each index i ∈ {1..n} to min(i, n+1−i). A symmetric array is fixed by its entries on the
octant cells (a, b) with 1 ≤ a ≤ b ≤ h; there are h(h+1)/2 = A008805(n−1) of them, and every choice of octant
entries gives exactly one symmetric array. So a(n) ≤ 2^(h(h+1)/2).

Call two octant cells *neighbors* if some cell of the array in the orbit of the first is a king-move neighbor of
some cell in the orbit of the second. Because every symmetry of the square maps king-move neighbors to king-move
neighbors, all cells of one orbit see the same multiset of orbits around them. So the condition "this 1 has a
neighboring 1" is the same for every cell of an orbit, and it holds iff some neighbor (in the octant sense, which
can be the cell itself, e.g. the middle cell when n is even) is 1. Neighbors are computed from the representative
cell (a, b) of each orbit. Cells outside the array are 0 and never help.

## 3. Transfer matrix (`tools/pinhole_tm.cpp`)

The octant cells are placed one at a time, in order of b, then a. Every neighbor of a cell lies within one octant
row of it, so a cell can be forgotten once its last neighbor has been placed. The *frontier* is the set of placed
cells that are not yet forgotten; for each, the state keeps 2 bits: its value, and (for a 1) whether a neighboring 1
has been seen yet. Placing a cell with value 1 marks it and all its placed 1-neighbors as satisfied. When a cell is
forgotten it must be 0 or satisfied, otherwise the partial array is dropped. The frontier has h + 1 cells, so a
state fits in one 64-bit word; states are kept in a hash table with their counts, and the counts are added.

The number of distinct states grows about 2.24-fold per unit of h (44,144,204 states at h = 22, 222,880,437 at
h = 24), far fewer than the 2^(h(h+1)/2) arrays.

Arithmetic: the counts are exact 384-bit integers up to n = 43 (enough, since a(n) < 2^300 for n ≤ 48). For
n = 43..47 a leaner build (`-DMODP`) keeps one 64-bit residue per state; each n is run modulo 2^64 and modulo the
four primes 2^62 − 57, 2^62 − 87, 2^62 − 117, 2^62 − 143 (`primes.txt`), and `crt.py` recombines them. The
moduli cover 312 bits and a(n) < 2^300 for n ≤ 48, so the reconstruction is exact.

| n | octant cells | max states | time | peak memory |
|---|---|---|---|---|
| 41 | 231 | 19,646,019 | 70 s (exact) | 11 GB |
| 43 | 253 | 44,144,204 | 273 s (exact); 51 s per residue | 22 GB; 6.7 GB |
| 45 | 276 | 99,191,126 | 2–3 min per residue | 13 GB |
| 47 | 300 | 222,880,437 | about 11 min per residue | 24 GB |

n = 48 would need about twice the memory of n = 47.

## 4. Checks

* **Brute force** (`tools/count_pinhole.c`): every symmetric array, by enumerating the octant, for n ≤ 14
  (2^28 arrays at n = 13, 14). Agrees with the transfer matrix for all n ≤ 14.
* **Exact against residues**: at n = 43 the CRT of the five residue runs equals the exact 384-bit count.
* **PARI** (`tools/check.gp`, the program in the entry): agrees for n = 1..10 (3 s).
* **By hand**: a(1) = 1 (the 1 × 1 array [1] has no neighbor), a(2) = 2 (all 0 or all 1), a(3) = 6 (of the 8
  choices of corner, edge and center, "center only" and "corners only" leave a lonely 1).
* **Without symmetry** the same rule gives A202788: the brute force `count_pinhole.c 1 n` reproduces its terms
  1, 12, 398, 45519, 21649928 for n = 1..5, which checks the reading of the rule (frame = 0, king adjacency).

## 5. Reproduce

```
g++ -O2 -march=native -o pinhole_tm tools/pinhole_tm.cpp
g++ -O2 -march=native -DMODP -o pinhole_tm_mod tools/pinhole_tm.cpp
./pinhole_tm 21                       # n, a(n), max states, frontier size
./pinhole_tm_mod 45 0                 # a(45) mod 2^64; repeat with each prime in primes.txt
python tools/crt.py 45                # reads tools/tm_mod_45.txt
gcc -O2 -o count_pinhole tools/count_pinhole.c && ./count_pinhole 0 13    # brute force
gp -q tools/check.gp                  # PARI, n = 1..10
```

Tracing the plates needs the Met photographs (public domain; object IDs 662882–662940, image files listed in
`tools/met1603_ids.txt`): `python tools/trace.py page.jpg page.npy`, then `tools/panels.py`, `tools/build_panels.py`, and
the rule tests `tools/rare.py`, `tools/props.py`, `tools/loo.py`.

## 6. A second rule: no small enclosed regions (`tools/growth_rules.py`, `tools/small_open.py`)

The pinhole is the smallest case of a more general picture: a group of open meshes walled in by darned meshes.
Call a *region* a set of open meshes joined through sides or corners, with the frame darned. In all 19 panels
every region has **at least 5 meshes**. The bound is tight: the only 5-mesh region is the open cross at the center of
the page 23 panel. The others have 9 meshes or many more. In 10 panels all the open meshes form a single region.

For comparison, `tools/null_sampler.c` draws 20 random symmetric no-pinhole designs per panel, with the same size
and the same fraction of darned meshes. They break the 5-mesh rule on most panels. On the two densest panels (pages
23 and 24) none of the 20 passes, and on the other panels 25–90% pass. Taken together, the chance that 19 random
designs would all pass is below 10^-8. The reason is that Vinciolo's open meshes form one web around the motifs,
while random designs have between 2 and 91 regions (median per panel), mostly tiny.

Rules on how a design grows fail. We tested 24 versions: "every open (or every darned) mesh touches a mesh of its
color that is closer to the center (or to the frame)", with square, diamond or round distance and side or corner
neighbors. Each version is broken in every one of the 19 panels, because the motifs (scrolls, hooks, flower heads)
curl back toward the center. On the darned side, islands of 2 darned meshes occur in two panels, so there is no
darned-region rule either.

The corner contacts matter: counted through sides only, the panels contain 1,792 open regions of fewer than 5
meshes, many of them single open meshes joined to the rest only at a corner, a standard stepped-diagonal device of
Vinciolo's.

## 7. The minimum-region table (`tools/minregion_tm2.cpp`)

Brute force (`tools/count_minregion.c`) counts T(n,k) for n ≤ 12. Beyond that, a transfer matrix runs over the same
octant cells as section 3, in the same order, with the same frontier. Two octant cells are joined if some cells of
their orbits are king-move neighbors. For each frontier cell the state stores which component it belongs to (0 for a
0-cell). For each component it stores:

* S, the sum of the orbit sizes (1, 4 or 8) of its octant cells;
* whether it touches the diagonal mirror line and whether it touches the vertical midline. A cell touches a line if
  its mirror image is itself or a king-move neighbor.

**Region sizes.** A region X of the full array that meets the octant contains every octant cell of its component C.
This is because crossing a mirror line by a king move makes X symmetric about that line. The symmetries that fix X
form the group H generated by the reflections in the lines C touches. Its order is 1, 2, 2 or 8 (none, the diagonal
only, the midline only, both). So

  |X| = S · |H| / 8,

and each orbit splits evenly among the 8/|H| copies of X. Once S · |H| ≥ 8k the region is large enough, so S is
capped there. When a component leaves the frontier for good, its size is checked against k.

A state packs 4-bit labels for up to 20 frontier cells and 11 bits (9 for the capped S, 2 for the flags) for up to
10 components into 192 bits. Counts are kept modulo 2^128 and modulo 2^61 − 1 and combined by CRT, which is exact
below 2^189 (the octant has at most 189 cells for n ≤ 37). `tools/minregion_tm.cpp` is a simpler first version with
string keys. `tools/run_table.py` runs the jobs and `tools/assemble_table.py` builds the b-file. Columns 1 and 2 are
not run: T(n,1) = 2^A008805(n−1), and column 2 is `b_pinhole_d4.txt`.

The hardest entries in each antidiagonal sit in its middle. The last two run, T(29,5) and T(27,7), took 7–8
minutes each with about 50 million states and up to 14 GB. Antidiagonal 35 needs larger jobs.

## 8. Checks for the table

* **Brute force** (`tools/count_minregion.c`, `tools/minregion_counts.txt`): all 96 entries with n ≤ 12, k ≤ 8.
* **Two transfer matrices**: `minregion_tm.cpp` and `minregion_tm2.cpp` agree on test entries up to n = 24.
* **Column 2** equals `b_pinhole_d4.txt` through n = 32. These terms pass 2^128, so the CRT path is tested.
* **PARI** (`tools/check_table.gp`, the program in the entry): agrees on antidiagonals 2..10 (45 terms).
* **Structure** (`tools/assemble_table.py`): every row is nonincreasing in k, and T(n,k) = 1 for k > n².

```
g++ -O2 -std=c++17 -o minregion_tm2 tools/minregion_tm2.cpp
./minregion_tm2 20 5                       # n k T(n,k)
python3 tools/run_table.py 34 4 4          # max n+k, workers, GB per job (resumable)
python tools/assemble_table.py             # writes b_minregion_table.txt
```

## Files

* `b_pinhole_d4.txt`: a(1)–a(47).
* `tools/pinhole_tm.cpp`, `tools/crt.py`, `tools/primes.txt`, `tools/tm_*.txt`: the transfer matrix and its runs.
* `tools/count_pinhole.c`, `tools/brute_13_14.txt`, `tools/check.gp`: the brute-force and PARI checks.
* `panels.txt`: the 19 traced panels. `tools/trace.py`, `panels.py`, `build_panels.py`, `export_panels.py`,
  `rules.py`, `rare.py`, `props.py`, `loo.py`: tracing and rule tests. `tools/met1603_ids.txt`: the Met objects and image files used.
* `b_minregion_table.txt`: T(n,k), antidiagonals 2..34. `tools/minregion_tm2.cpp`, `tools/minregion_tm.cpp`,
  `tools/run_table.py`, `tools/assemble_table.py`, `tools/table_runs.txt`: the table's transfer matrix and runs.
  `tools/count_minregion.c`, `tools/minregion_counts.txt`, `tools/check_table.gp`: its checks.
* `tools/growth_rules.py`, `tools/small_open.py`, `tools/null_sampler.c`, `tools/growth_load.py`: the region and
  growth rule tests and the random comparison designs.
* `NOTES.md`: the working log.
