# A119501 in plain English

*How many polyhedra have n corners, if a shape and its mirror image count as two?*

---

## 1. The idea

This is the [main page's](../layman.md) question, how many polyhedron shapes have n corners, with one
change. There, a shape and its mirror image are the same shape. Here, you may turn a polyhedron around
but not reflect it, so a *chiral* shape (one that differs from its mirror image, like a left and a right
hand) counts twice: once for each version.

## 2. Where this shows up

* Chemistry: left- and right-handed versions of a molecule's shape can behave differently, so a
  catalogue of shapes "up to rotation" is the relevant one.
* Counting theory: this "turning only" count is the natural first step; the main sequence A000944 is
  built from it plus a count of mirror-symmetric shapes.

## 3. What A119501 counts

The number of polyhedra with n corners, mirror images counted separately.

* 4 corners: the tetrahedron, 1. 5 corners: the square pyramid and the triangular double pyramid, 2.
  All of these are their own mirror images.
* 6 corners: there are 7 polyhedron shapes, and exactly one is chiral: a triangular prism with one of
  its square faces cut into two triangles along a diagonal. The diagonal can slope "left" or "right",
  and no turning changes one into the other. So this count is 7 + 1 = 8.

      n:      4  5  6   7    8     9      10 ...              17
      a(n):   1  2  8  45  419  4798   62754 ... 12831576165782

(The entry starts at n = 1; the terms for 1 to 3 corners are 0.)

## 4. Why it's hard

The numbers grow about 17-fold per extra corner, like the main sequence. The known values stopped at
17 corners, about 10^13 shapes; methods that build and inspect the shapes one at a time cannot go much
further.

## 5. What we found

Twenty-three new terms, from 18 to 40 corners:

**a(18) = 215,707,912,782,073**, **a(19) = 3,675,708,943,760,861**, …,
**a(40) = 1,360,698,139,575,627,018,329,534,273,611,152,495,438,033** (about 1.4·10^42).

They come from the [A000944](../layman.md) computation, which counts polyhedra by corners and faces with
formulas and Burnside's lemma, separately for "turning only" and "turning or reflecting". This sequence
is the "turning only" count added up over all numbers of faces. That half of the computation needs only
rotations, so it is the easier half.

For large n almost every polyhedron is chiral, so this count is very nearly twice A000944. At 40
corners, twice A000944(40) minus a(40) leaves 9,122,358,895,023,892,499,374,815: the number of
mirror-symmetric polyhedra with 40 corners, a tiny fraction of the 6.8·10^41 shapes.

## 6. How we know it's right

* All 17 previously known terms come out exactly.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; every rotation
  count from the formulas matched.
* The same "turning only" table, read along a different line, gives the count of triangle-only
  polyhedra with mirror images separate (A253882, from an independent older formula), and it matches
  for every n up to 40.
* Swapping corners and faces must leave the full table unchanged, and it does.

## 7. Why it matters (honestly)

This is an extension of a standard companion table, not a breakthrough.
* It goes from 17 to 40 corners, matching the range of the main sequence.
* It is the easier half of the main computation; its value is completeness, so that the "mirror images
  separate" and "mirror images identified" counts by corners now reach equally far.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Chiral:** different from its mirror image, like a hand.
- **Up to rotation:** two shapes count as the same if one can be turned to match the other, but not
  reflected.
