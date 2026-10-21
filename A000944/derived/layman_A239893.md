# A239893 in plain English

*The table of polyhedra by faces and corners, with left- and right-handed versions counted apart.*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

Think of a pair of gloves: mirror images that no amount of turning makes identical. Some polyhedra are
like that. A239893 counts shapes with n faces and k corners, where a shape and its mirror image count
as two unless you can *turn* one into the other. It is the "handed" version of the table
[A212438](layman_A212438.md), with the same layout.

## 2. A "triangle read by rows"

The OEIS stores a table as one long list, reading the rows one after another. Row n lists the shapes
with n faces and k = 4, 5, ..., 2n − 4 corners, with 0 where none exist. Next to each row, the same row
when mirror images are *not* counted apart (A212438):

    faces n = 4:  1                      (1)
    faces n = 5:  0  1  1                (0  1  1)
    faces n = 6:  0  1  3  2  2          (0  1  2  2  2)
    faces n = 7:  0  0  2 11 16 10  6    (0  0  2  8 11  8  5)

Where the numbers differ, there are shapes with no mirror symmetry. The first one is row 6, k = 6:
of the two shapes with 6 faces and 6 corners, the pentagonal pyramid is its own mirror image, but the
other comes in a left and a right version. So the handed count is 1 + 2 = 3.

In general the handed number is somewhere between the plain number and twice it. For large tables it
is almost exactly twice, because almost every big polyhedron has no symmetry at all.

## 3. Swapping corners and faces

Putting a corner in the middle of each face and joining neighbours gives the *dual* polyhedron (cube ↔
octahedron), and it keeps handedness. So the number with n faces and k corners equals the number with
k faces and n corners: row 6, k = 7 and row 7, k = 6 are both 2, and row 7, k = 8 and row 8, k = 7 are both 16.

## 4. What we found

The entry's data file covered rows 4 to 17 (196 numbers). We extend it to 40 faces: **1369 numbers**,
rows 18 to 40 new. For example:

* 18 faces and 20 corners: **10,541,534,121,787** (against 5,270,785,332,349 without handedness; just
  under twice as many)
* 40 faces and 76 corners: **286,312,125,410,688,368,031,270,157,720**

The numbers come from the computation on the main page. There the handed count is the easier half: it
needs only the rotation symmetries, and the mirror symmetries are added afterwards for A212438.

## 5. How we know it's right

* All 196 previously known numbers (rows 4 to 17, computed by Andrew Howroyd) come out exactly.
* Every polyhedron with up to 11 corners was built by computer with its full list of symmetries, and
  the formula counts matched them.
* The table is symmetric under swapping corners and faces, as it must be.
* Combined with the mirror counts, the same numbers give the plain table A212438 and the main polyhedra
  count, which match all their own known values.

## 6. Why it matters (honestly)

This is a small result: a standard refinement table, now much longer.
* The value is completeness: 1369 exact entries instead of 196, out to 40 faces.
* Its middle diagonal is [A342057](layman_A342057.md), and its row sums count polyhedra up to turning
  (A119501); both extend with it.
* Together with A212438 it tells how many shapes of each size are mirror-symmetric.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Chiral:** different from its mirror image, like a left glove.
- **Sensed:** the technical word for "counted up to turning only", so mirror images are distinct.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
