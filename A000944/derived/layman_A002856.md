# A002856 in plain English

*How many polyhedra have exactly as many faces as corners?*

---

## 1. The idea

A polyhedron is a solid with flat faces: a cube, a pyramid, a soccer ball. As on the
[main page for polyhedra](../layman.md), two polyhedra count as the same *shape* if they are built the
same way (same corners, joined by edges in the same pattern), whatever the lengths and angles.

Most polyhedra have different numbers of corners and faces. A cube has 8 corners but only 6 faces. Some
are balanced: they have exactly as many faces as corners. This sequence counts the balanced ones.

## 2. Where this shows up

* Every pyramid is balanced. A pyramid over a polygon with m sides has m + 1 corners (the base plus the
  tip) and m + 1 faces (m triangles plus the base).
* Swapping corners and faces (put a corner in the middle of each face, join neighbours) turns every
  polyhedron into another one, its *dual*. The cube becomes the octahedron. Balanced polyhedra are
  exactly the ones whose dual is again balanced, so they sit in the middle of that swap.
* This is the middle diagonal of the big table of polyhedra by corners and faces
  ([A212438](layman_A212438.md)).

## 3. What A002856 counts

The number of polyhedron shapes with n corners and n faces.

* n = 4: only the tetrahedron (4 corners, 4 triangular faces).
* n = 5: only the square pyramid (5 corners: four on the base and the tip; 5 faces: four triangles and
  the square). It has 8 edges, and it is the only shape with 5 and 5.
* n = 6: two shapes, one of them the pentagonal pyramid.

        n:      4  5  6  7   8    9    10     11 ...            18
        a(n):   1  1  2  8  42  296  2635  25626 ... 700335433295

## 4. Why it's hard

The counts grow 12- to 14-fold for each extra corner. The known values (up to 18) came from programs
that build every polyhedron once and sort them by corners and faces. Past 18 corners there are too many
polyhedra to build one by one.

## 5. What we found

All terms from n = 19 to n = 40. For example:

* **a(19) = 8,816,236,975,026**
* **a(20) = 112,504,197,045,246**
* **a(40) = 7,008,894,109,272,003,406,002,166,618,568,677,332** (about 7·10^36)

They come out of the same computation as the main page: polyhedra are counted by formulas (with
Burnside's lemma handling the symmetric ones), not built one at a time, and the counts are kept separate
for every number of corners and faces. This sequence is one diagonal line of that table.

## 6. How we know it's right

* All 15 previously known terms (n = 4 to 18) come out exactly. These were found by building the shapes,
  a completely different method.
* For up to 11 corners we built every polyhedron by computer and checked the formula counts against
  them, symmetry type by symmetry type.
* The whole table is symmetric under swapping corners and faces, as it must be, even though the method
  treats corners and faces differently.

## 7. Why it matters (honestly)

This is a small result. A002856 is a standard "refinement" of the count of all polyhedra.
* It more than doubles the entry, from 18 to 40, and the new terms are exact.
* The value is completeness: the balanced count now goes as far as the main polyhedra count.
* No new idea was needed here; it is a by-product of the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
- **Burnside's lemma:** a rule for counting things "up to symmetry" by averaging over the symmetries.
