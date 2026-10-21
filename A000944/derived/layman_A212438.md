# A212438 in plain English

*The table of polyhedra by faces and corners, with the impossible cells written as zeros.*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

A212438 answers: how many shapes have n faces and k corners? It holds the same numbers as
[A058787](layman_A058787.md), laid out differently. Every row runs over the same simple range of k, from
4 corners up to 2n − 4 corners, and writes 0 where no polyhedron exists.

## 2. A "triangle read by rows"

The OEIS stores a table as one long list, reading the rows one after another. Row n lists the
polyhedra with n faces and k = 4, 5, ..., 2n − 4 corners:

    faces n = 4:  k = 4                    1
    faces n = 5:  k = 4, 5, 6              0  1  1
    faces n = 6:  k = 4 ... 8              0  1  2  2  2
    faces n = 7:  k = 4 ... 10             0  0  2  8 11  8  5

So the list starts 1, 0, 1, 1, 0, 1, 2, 2, 2, 0, 0, 2, 8, ... Things to check by hand:

* Row 5, k = 4 is 0: only the tetrahedron has 4 corners, and it has 4 faces, not 5.
* Row 5, k = 5 and k = 6: the square pyramid and the triangular prism.
* The zeros at the start of each row are shapes with too few corners for that many faces (each face
  needs at least 3 sides, which forces at least n/2 + 2 corners).

The rows stop at 2n − 4 because a polyhedron with n faces has at most that many corners (the cube,
6 faces and 8 corners, is an example at the limit; [A058787](layman_A058787.md) explains why).

## 3. Swapping corners and faces

Put a new corner in the middle of each face and join neighbours: you get the *dual* polyhedron (the cube
becomes the octahedron). So the number with n faces and k corners equals the number with k faces and
n corners: row 6, k = 7 and row 7, k = 6 are both 2. This symmetry is a free check on every number.

## 4. What we found

The entry's data file covered rows 4 to 17 (196 numbers). We extend it to 40 faces: **1369 numbers**,
rows 18 to 40 new. For example:

* 18 faces and 20 corners: **5,270,785,332,349** shapes
* 19 faces and 19 corners: **8,816,236,975,026**
* 40 faces and 76 corners: **143,156,062,705,641,154,797,382,062,151**

Row 40 adds up to the main page's number of polyhedra with 40 corners (by the swap above), about 7·10^41.
The numbers come from the table computed for the main page, by formulas rather than by building shapes.

## 5. How we know it's right

* All 196 previously known numbers (rows 4 to 17, computed by Andrew Howroyd) come out exactly.
* Every polyhedron with up to 11 corners was built by computer, and the formula counts matched them.
* The table is symmetric under swapping corners and faces (checked on 324 pairs), although the method
  treats corners and faces differently.
* The zeros land exactly where the corner limits say they must.

## 6. Why it matters (honestly)

This is a small result: a standard refinement table, now much longer.
* The value is completeness: 1369 exact entries instead of 196, out to 40 faces.
* Several other OEIS sequences are lines through this table (balanced polyhedra
  [A002856](layman_A002856.md), [A058789](layman_A058789.md), triangulations), so they all extend at once.
* No new idea was needed; it is a by-product of the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Triangle read by rows:** a table with rows of growing length, stored as one list row after row.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
