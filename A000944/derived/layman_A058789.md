# A058789 in plain English

*How many polyhedra have one more corner than faces?*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

This sequence asks for polyhedra that are *almost* balanced: n faces and n + 1 corners. Its neighbour
[A002856](layman_A002856.md) counts the exactly balanced ones (n faces, n corners).

## 2. Where this shows up

* Swapping corners and faces turns every polyhedron into another one, its *dual* (cube ↔ octahedron).
  So "n faces and n + 1 corners" and "n corners and n + 1 faces" give the same count: each shape of one
  kind is paired with exactly one of the other.
* In the table of polyhedra by corners and faces ([A212438](layman_A212438.md)) this is the line just
  next to the middle diagonal.

## 3. What A058789 counts

The number of polyhedron shapes with n faces and n + 1 corners.

* n = 4: none. The only polyhedron with 4 faces is the tetrahedron, which has 4 corners, not 5.
* n = 5: one, the triangular prism (a Toblerone box): 2 triangles and 3 rectangles make 5 faces, and it
  has 6 corners. Its dual, the triangular double pyramid, has 5 corners and 6 faces.
* n = 6: two shapes with 6 faces and 7 corners.

        n:      4  5  6   7   8    9    10     11 ...             18
        a(n):   0  1  2  11  74  633  6134  64439 ... 2173270387051

## 4. Why it's hard

The counts grow 12- to 15-fold at each step. The known values (up to 18 faces) came from programs that
build every polyhedron once; beyond that there are too many to build.

## 5. What we found

All terms from n = 19 to n = 40. For example:

* **a(19) = 27,777,985,734,886**
* **a(20) = 359,289,661,461,953**
* **a(40) = 25,220,556,563,347,336,413,383,575,055,259,960,478** (about 2.5·10^37)

They come from the same table as the main page: polyhedra counted by formulas for every number of
corners and faces, with Burnside's lemma handling the symmetric ones, and no shape built one at a time.
This sequence is one line of that table.

## 6. How we know it's right

* All 15 previously known terms (n = 4 to 18) come out exactly. Those were found by building every
  shape, a completely different method.
* For up to 11 corners we built every polyhedron by computer and checked the formula counts against
  them, symmetry type by symmetry type.
* "n faces, n + 1 corners" and "n corners, n + 1 faces" are computed separately and agree every time,
  as duality demands, even though our method treats corners and faces differently.

## 7. Why it matters (honestly)

This is a small result: a standard refinement of the count of all polyhedra.
* The entry goes from 18 to 40 terms, all exact.
* Its value is completeness, matching the range of the main polyhedra count.
* No new idea was needed; it is a by-product of the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
- **Burnside's lemma:** a rule for counting things "up to symmetry" by averaging over the symmetries.
