# A000109 in plain English

*How many different polyhedra can be built entirely out of triangles?*

---

## 1. The idea

Take a polyhedron whose faces are all triangles: the tetrahedron, the octahedron, the icosahedron, or a
lumpy shape like a geodesic dome. Mathematicians call these *simplicial* polyhedra. As on the
[main page](../layman.md), only the pattern matters: which corners are joined by edges. Lengths and
angles are forgotten, so a squashed octahedron is still the octahedron.

Another way to see the same objects: draw dots on a sphere and join them by curves into triangles that
cover the whole sphere. That is a *triangulation of the sphere*, and it is the same thing.

## 2. Where this shows up

* Computer graphics and engineering: surfaces are stored as triangle meshes; these are the closed ones.
* Chemistry: the "deltahedral" cages of boron clusters (boranes) are triangulated spheres.
* Graph theory: these are exactly the *maximal planar graphs*: networks drawn in the plane without
  crossings, to which no further edge can be added.
* Mirror partners: swap corners and faces (put a dot in each triangle, join dots of neighbouring
  triangles) and every triangle polyhedron becomes a polyhedron where every corner meets exactly 3 edges.
  The octahedron (6 corners) turns into the cube (8 corners). So A000109(n) also counts these "cubic"
  polyhedra with 2n − 4 corners, which is why the entry's name mentions both.

## 3. What A000109 counts

The number of triangle-only polyhedra with n corners.

* 4 corners: only the tetrahedron.
* 5 corners: only the triangular double pyramid (two tetrahedra glued face to face).
* 6 corners: two. The octahedron, where every corner meets 4 edges; and the double pyramid with a
  tetrahedron glued onto one face, where the corners meet 3, 3, 4, 4, 5 and 5 edges. Different counts
  of edges at the corners means different shapes.

Swap corners and faces and these become the tetrahedron, the triangular prism (6 corners), and the cube
plus one more 8-cornered shape: 2n − 4 corners each time.

    n:      3  4  5  6  7   8   9    10    11 ...              24
    a(n):   1  1  1  2  5  14  50   233  1249 ... 230214772514389

(The first term, n = 3, counts a single triangle seen as two faces glued at the rim, by convention.)

## 4. Why it's hard

The numbers grow about 8-fold per extra corner. The known values come from Brinkmann and McKay's program
plantri, which builds every triangulation once. It reached 24 corners, about 2.3·10^14 shapes; going
further one at a time is not practical.

There is a twist. If you count shapes *as seen in a mirror as different* (a left-handed and a
right-handed version count twice), formulas exist, and that version (A253882) was already known to 256
terms. What was missing is the count where a shape and its mirror image are the same. That needs the
number of shapes which have a mirror symmetry, and that was the hard part.

## 5. What we found

Fifty-six new terms, from 25 to 80 corners:

**a(25) = 1,865,788,057,030,046**, **a(26) = 15,223,573,236,152,171**, …,
**a(40) = 143,156,062,705,641,154,797,382,062,151** (about 1.4·10^29), …,
**a(80) = 13,610,262,362,737,077,739,549,176,947,869,565,030,337,043,689,010,224,527,879,637,230,383** (about 1.4·10^67).

These came for free from the computation behind [A000944](../layman.md): we counted polyhedra by corners
*and* faces, and a triangle polyhedron with n corners is exactly one with 2n − 4 faces. The new ingredient
is the one described there: counting mirror-symmetric shapes by cutting them along the mirror and
counting the halves. That is exactly the piece this sequence lacked.

Terms 41 to 80 use a shortcut. Every triangle polyhedron is counted once if it is its own mirror image
and twice (left and right versions) otherwise. So the plain count is the mirror-images-different count
(A253882, from its formula) plus the number of mirror-symmetric shapes, all divided by two. For
triangle-only shapes the mirror halves are simple: each is a triangulated patch whose rim records
where the mirror passes through a corner or cuts across an edge. Counting those patches takes seconds
instead of hours. The rarer "turn-and-reflect" symmetries still use the full method.

## 6. How we know it's right

* All 22 previously known terms (n = 3 to 24, from plantri) come out exactly.
* The same computation, done without mirrors, gives A253882 (the count where mirror images differ), and
  it matches that entry for every n up to 40. A253882 comes from an older, independent formula.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; every symmetry
  count from our formulas matched. This rules out a wrong formula for any kind of symmetry.
* Swapping corners and faces must give the same totals, and it does, even though our mirror method
  treats corners and faces differently.
* The shortcut used for 41 to 80 reproduces every term up to 40 exactly.
* The numbers grow as theory predicts: each new term is about 256/27 ≈ 9.48 times the last, with a
  correction that fades as n grows, and our ratios approach that value smoothly (9.435 at n = 41,
  9.470 at n = 80).

## 7. Why it matters (honestly)

This is an extension of a classic sequence, not a new theorem.
* A000109 is one of the OEIS's "core" sequences; the new terms are the first beyond what exhaustive
  generation reached, and they take it from 24 to 80 corners.
* It shows the gap clearly: counting with mirrors-as-different was easy and done long ago; counting with
  mirrors-as-same is what the new mirror method makes possible.
* The method is not limited by how many shapes exist, so further terms only cost more computer time.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Simplicial polyhedron / triangulation of the sphere:** a polyhedron all of whose faces are triangles.
- **Cubic polyhedron:** one in which every corner meets exactly 3 edges; the "corners ↔ faces" partner of
  a triangle polyhedron.
- **Mirror image (chirality):** a shape is *chiral* if its mirror image is a different shape, like a left
  and a right hand.
