# A342057 in plain English

*Balanced polyhedra (as many faces as corners), with left- and right-handed versions counted apart.*

---

## 1. The idea

A polyhedron is a solid with flat faces. As on the [main page for polyhedra](../layman.md), two
polyhedra count as the same *shape* if their corners are joined by edges in the same pattern, whatever
the lengths and angles.

Think of a pair of gloves. The left glove and the right glove are mirror images: same pattern, but no
amount of turning makes one into the other. Some polyhedra are like that too. Most counts (such as
[A002856](layman_A002856.md)) treat a shape and its mirror image as one shape. This one counts them
separately: two shapes are the same only if you can *turn* one into the other, without a mirror.

A342057 does this for balanced polyhedra: n corners and n faces.

## 2. Where this shows up

* Chemistry: molecules that differ only by being mirror images ("chiral" molecules) can behave very
  differently, so counting with handedness matters there.
* Comparing the two counts tells how many shapes have a mirror symmetry: twice A002856 minus A342057
  is exactly the number of balanced shapes that are their own mirror image.
* This is the middle diagonal of the handed table of polyhedra by corners and faces
  ([A239893](layman_A239893.md)).

## 3. What A342057 counts

The number of polyhedron shapes with n corners and n faces, where a shape and its mirror image count
as two if no rotation turns one into the other.

* n = 4: the tetrahedron, 1. n = 5: the square pyramid, 1. Both are their own mirror images.
* n = 6: there are two shapes with 6 corners and 6 faces. One is the pentagonal pyramid, which is its
  own mirror image. The other has no mirror symmetry, so it comes in a left and a right version. Total:
  1 + 2 = **3**, against 2 in A002856.

        n:      4  5  6   7   8    9    10     11 ...             17
        a(n):   1  1  3  11  69  541  5096  50586 ... 112984297173

## 4. Why it's hard

The counts grow 12- to 14-fold at each step. The known values (up to 17) came from building and
sorting every polyhedron; beyond that there are too many.

## 5. What we found

All terms from n = 18 to n = 40. For example:

* **a(18) = 1,400,659,574,748**
* **a(19) = 17,632,426,037,198**
* **a(40) = 14,017,788,218,544,004,385,875,816,658,241,699,960** (about 1.4·10^37)

For large n this is almost exactly twice A002856, because almost every large polyhedron has no symmetry
at all, so it comes in two mirror versions. For n = 40 only
2,426,128,516,578,895,654,704 of the balanced shapes are their own mirror image, about one in 3·10^15.

The method is the one on the main page. Counting "up to turning only" is actually the easier half of
that computation: it needs only the rotation symmetries. The mirror symmetries are needed for A002856.

## 6. How we know it's right

* All 14 previously known terms (n = 4 to 17) come out exactly.
* For up to 11 corners we built every polyhedron by computer, listed its symmetries, and checked every
  formula count against them.
* The handed table is symmetric under swapping corners and faces, as it must be.
* Combined with the mirror counts, the same numbers give A002856 and the main polyhedra count, which
  match their own known terms.

## 7. Why it matters (honestly)

This is a small result: a standard refinement table, read along its diagonal.
* The entry goes from 17 to 40, all exact; its value is completeness.
* Paired with A002856 it gives the number of mirror-symmetric balanced polyhedra up to 40 corners.
* No new idea was needed; it is a by-product of the main computation.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Chiral:** different from its mirror image, like a left glove.
- **Burnside's lemma:** a rule for counting things "up to symmetry" by averaging over the symmetries.
