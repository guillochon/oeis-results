# A058787 in plain English

*A table of polyhedra sorted by faces and corners.*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

The main page counts shapes by corners alone. A058787 sorts them more finely: how many shapes have
n faces *and* k corners? The answers form a table, one row for each number of faces.

## 2. A "triangle read by rows"

The OEIS stores tables like this as one long list, reading the rows one after another. Here row n lists
the polyhedra with n faces, for every possible number of corners k from smallest to largest:

    faces n = 4:  k = 4                     1
    faces n = 5:  k = 5, 6                  1  1
    faces n = 6:  k = 5, 6, 7, 8            1  2  2  2
    faces n = 7:  k = 6, 7, 8, 9, 10        2  8 11  8  5

So the sequence starts 1, 1, 1, 1, 2, 2, 2, 2, 8, 11, 8, 5, ... The rows get longer, so the picture looks
like a triangle. Some entries to check by hand:

* Row 4: only the tetrahedron has 4 faces.
* Row 5: the square pyramid (5 corners) and the triangular prism (6 corners).
* Row 6, k = 8: the cube and one other shape with 6 faces and 8 corners.

Adding up a row gives the number of polyhedra with n faces, which is the same as the main page's count
of polyhedra with n corners (see below).

## 3. Why the rows stop where they do

With n faces, a polyhedron has at most 2n − 4 corners. In plain words: every corner has at least 3 edges
coming out of it, and every edge has 2 ends, so there are at least one and a half times as many edges
as corners. Euler's rule says corners minus edges plus faces is always 2. Put the two together and the
corners can be at most twice the faces, minus 4. The cube is an extreme case: 6 faces, 8 corners, three
edges at every corner. In the same way, because every face has at least 3 sides, there are at least
n/2 + 2 corners.

The entry's title says the rows run up to k = 2n + 8. That is a typo (with 6 faces there can never be
20 corners); the rows themselves already stop at 2n − 4. Our proposed edit corrects the title.

## 4. Swapping corners and faces

Put a new corner in the middle of each face and join neighbours: you get another polyhedron, the
*dual*. The cube becomes the octahedron. Corners and faces trade places, so the number of shapes with
n faces and k corners equals the number with k faces and n corners. In the table: row 6, k = 7 is 2,
and so is row 7, k = 6. Row 6, k = 8 is 2 (cube and friend), and so is row 8, k = 6 (octahedron and
friend). That is why row sums equal the count by corners.

## 5. What we found

The entry listed the first 54 numbers (rows 4 to 11 and part of row 12). We give all 37 rows from 4 faces to 40 faces:
**1027 numbers**. Rows up to 17 faces also appear in the layout of [A212438](layman_A212438.md); rows 18
to 40 are new. For example:

* 18 faces and 20 corners: **5,270,785,332,349** shapes
* 40 faces and 40 corners: **7,008,894,109,272,003,406,002,166,618,568,677,332**
* 40 faces and 76 corners (the most possible): **143,156,062,705,641,154,797,382,062,151**

These all come from the table computed for the main page, by formulas rather than by building shapes.

## 6. How we know it's right

* All 54 previously listed numbers come out exactly, and so do all 196 numbers of the same table in
  A212438 (up to 17 faces), which were found by a different method.
* Every polyhedron with up to 11 corners was built by computer, and the formula counts matched.
* The table is symmetric under swapping corners and faces (checked on 324 pairs), although the method
  treats corners and faces differently.
* Each row adds up to the known count of all polyhedra (up to 18).

## 7. Why it matters (honestly)

This is a small result: a standard refinement table, now much longer.
* The value is completeness: 1027 exact entries instead of 54, out to 40 faces.
* It fixes a typo in the entry's title.
* No new idea was needed; it is a by-product of the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Euler's rule:** for any polyhedron, corners − edges + faces = 2.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
