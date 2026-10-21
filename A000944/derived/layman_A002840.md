# A002840 in plain English

*How many different polyhedra have exactly n edges?*

---

## 1. The idea

A polyhedron is a solid with flat faces, like a cube or a pyramid. As on the
[main page](../layman.md), two polyhedra count as the same shape if their corners, edges and faces are
connected in the same pattern, whatever the lengths and angles. The main page sorts shapes by number of
corners; this sequence sorts them by number of *edges* instead.

## 2. Where this shows up

* Graph theory: by Steinitz's theorem these shapes are exactly the "3-connected planar graphs", so this
  counts those networks by number of links.
* It is the standard companion table to the counts by corners (A000944) and by faces; catalogues of
  small polyhedra are usually organised this way.

## 3. What A002840 counts

The number of polyhedron shapes with n edges.

* 6 edges: the tetrahedron.
* 7 edges: none. No polyhedron has exactly 7 edges. (Each face has at least 3 edges and each edge
  borders 2 faces, so 7 edges allow at most 4 faces; likewise at most 4 corners. But Euler's rule
  says corners + faces = edges + 2 = 9, and 4 + 4 is only 8.)
* 8 edges: the square pyramid.
* 9 edges: two, the triangular prism and the triangular double pyramid.
* 10 edges: two, the pentagonal pyramid and a triangular prism with one square face cut along a diagonal.

      n:      6  7  8  9  10  11  12  13  14 ...            35
      a(n):   1  0  1  2   2   4  12  22  58 ... 5986979643542

## 4. Why it's hard

The numbers grow by a factor of about 3.6 per extra edge. The known values came from computer programs
that build every polyhedron once and count them, and they reached 35 edges: about 6·10^12 shapes. Going much further one at a time is not practical.

## 5. What we found

Forty-four new terms, from 36 to 79 edges:

**a(36) = 21,636,767,535,118**, **a(40) = 3,792,268,711,213,413**, …,
**a(79) = 101,256,857,385,646,691,262,508,879,784,636,540,584** (about 10^38).

No new idea was needed beyond the [A000944](../layman.md) computation, which counts polyhedra by corners
and faces together using formulas and Burnside's lemma (no shapes are built). By Euler's rule a shape
with n edges has corners + faces = n + 2, so we add up the right entries of that table. Our table goes
up to 40 corners; with 79 edges, corners + faces = 81, so one of the two is at most 40, and swapping
corners with faces covers the other case. That is why we stop at 79.

## 6. How we know it's right

* All 30 previously known terms (6 to 35 edges, found by exhaustive generation) come out exactly.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; every count from
  the formulas matched.
* Swapping corners and faces must leave the counts unchanged, and it does, even though our method for
  mirror-symmetric shapes treats corners and faces differently.
* The same table reproduces many other published counts (by corners, by faces, by corners and faces
  together), so an error would have to hide from all of them at once.

## 7. Why it matters (honestly)

This is a straightforward extension of a standard table, not a discovery.
* It more than doubles the range of a reference sequence: from 35 to 79 edges.
* It also confirmed a correction: the companion sequence A005645 (mirror images counted separately)
  had one wrong digit at 34 edges, and this sequence's value at 34 is what pins it down
  (see [A005645](layman_A005645.md)).

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Euler's rule:** for any polyhedron, corners − edges + faces = 2.
- **Polyhedral graph:** the network of corners and edges of a polyhedron.
