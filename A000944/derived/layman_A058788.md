# A058788 in plain English

*A table of polyhedra sorted by edges and corners.*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

A058788 asks: how many shapes have n edges and k corners? Once you know the edges and corners, the
faces are fixed too, by Euler's rule (corners − edges + faces = 2). A cube has 8 corners and 12 edges,
so it must have 12 + 2 − 8 = 6 faces.

## 2. A "triangle read by rows"

The OEIS stores a table as one long list, reading the rows one after another. Row n lists the shapes
with n edges, for each possible number of corners k:

    edges n = 6:   k = 4              1
    edges n = 7:   (none)
    edges n = 8:   k = 5              1
    edges n = 9:   k = 5, 6           1  1
    edges n = 10:  k = 6              2
    edges n = 11:  k = 6, 7           2  2
    edges n = 12:  k = 6, 7, 8        2  8  2

So the list starts 1, 1, 1, 1, 2, 2, 2, 2, 8, 2, ... Things to check by hand:

* 6 edges: the tetrahedron. 8 edges: the square pyramid.
* 9 edges: the triangular double pyramid (5 corners) and the triangular prism (6 corners).
* 7 edges: no polyhedron at all. This is a famous little fact: every number of edges from 6 upward is
  possible except 7.
* The row ranges: every corner has at least 3 edges and every edge has 2 ends, so 3 × corners is at most
  2 × edges. That caps k at two thirds of n. The same argument for faces (each has at least 3 sides)
  gives at least n/3 + 2 corners.

## 3. Swapping corners and faces

Put a new corner in the middle of each face and join neighbours: you get the *dual* polyhedron, with the
same number of edges and corners and faces swapped (cube ↔ octahedron, both with 12 edges). So each
row reads the same forwards and backwards. In the 12-edge row, the first 2 (6 corners: the octahedron
and one other) mirrors the last 2 (6 faces: the cube and one other), and the 8 in the middle counts
shapes with 7 corners and 7 faces.

## 4. What we found

The entry listed the first 74 numbers (up to 26 edges). We give every row up to 79 edges:
**925 numbers**. For example, the row for 79 edges has 24 entries, from 29 to 52 corners:

* 79 edges, 29 corners: **3,466,081,847,146,905,007,622**
* 79 edges, 40 corners (the two middle entries, 40 and 41 corners):
  **25,220,556,563,347,336,413,383,575,055,259,960,478** each

That middle number is also the 40th term of [A058789](layman_A058789.md): 40 corners and 41 faces, or the
reverse, which is the same thing by the swap above. The row sums count all polyhedra by edges (A002840).

Why 79? The underlying table was computed for up to 40 corners (and 76 faces). Every shape with at most
79 edges has 40 corners or fewer, or 40 faces or fewer, so by the swap all of its rows are covered.

## 5. How we know it's right

* All 74 previously listed numbers come out exactly. So do the published counts of polyhedra by edges
  (A002840, up to 35 edges), which were found by building every shape.
* Every polyhedron with up to 11 corners was built by computer, and the formula counts matched them.
* Every row is a palindrome, as the corner-face swap demands, even though the method treats corners
  and faces differently.

## 6. Why it matters (honestly)

This is a small result: a standard refinement table, now much longer.
* The value is completeness: 925 exact entries instead of 74, up to 79 edges.
* It is the same information as the faces-and-corners table ([A058787](layman_A058787.md)),
  rearranged; no new idea was needed beyond the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Euler's rule:** for any polyhedron, corners − edges + faces = 2.
- **Dual polyhedron:** the shape you get by swapping the roles of corners and faces.
