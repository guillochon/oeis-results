# Symmetric lace designs without pinholes

This folder computes a(1)–a(47) of a new OEIS sequence:

> **a(n)** = number of n × n binary arrays that are invariant under the 8 symmetries of the square, in which every
> 1 is adjacent to another 1 horizontally, vertically or diagonally.

1, 2, 6, 7, 48, 55, 705, 861, 21454, 26696, 1311647, 1650868, 159700301, 203876423, … (`b_pinhole_d4.txt` has
n = 1..47). Without the symmetry, the same condition gives A202788.

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

## Files

* `b_pinhole_d4.txt`: a(1)–a(47).
* `tools/pinhole_tm.cpp`, `tools/crt.py`, `tools/primes.txt`, `tools/tm_*.txt`: the transfer matrix and its runs.
* `tools/count_pinhole.c`, `tools/brute_13_14.txt`, `tools/check.gp`: the brute-force and PARI checks.
* `panels.txt`: the 19 traced panels. `tools/trace.py`, `panels.py`, `build_panels.py`, `export_panels.py`,
  `rules.py`, `rare.py`, `props.py`, `loo.py`: tracing and rule tests. `tools/met1603_ids.txt`: the Met objects and image files used.
* `NOTES.md`: the working log.
