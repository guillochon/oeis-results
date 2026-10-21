# A006445 in plain English

*How many polyhedra with n edges look the same in a mirror?*

---

## 1. The idea

Hold a cube up to a mirror: the reflection is a cube. Hold up your right hand: the reflection is a left
hand, which no turning will make into a right hand. Shapes of the first kind are *achiral* (they have a
mirror symmetry); shapes of the second kind are *chiral*.

As on the [main page](../layman.md), a polyhedron's shape is the pattern of its corners, edges and faces,
forgetting lengths and angles. This sequence counts only the achiral shapes, sorted by number of edges.

## 2. Where this shows up

* Chemistry: whether a molecule's shape equals its mirror image decides whether it can come in left- and
  right-handed forms.
* Counting theory: it is the bridge between the two other edge counts. If S counts shapes with mirror
  images separate ([A005645](layman_A005645.md)) and U counts shapes with mirror images identified
  ([A002840](layman_A002840.md)), then this sequence is 2U − S.

## 3. What A006445 counts

The number of mirror-symmetric polyhedra with n edges.

* 6 edges: the tetrahedron. 8 edges: the square pyramid. Both are their own mirror images.
* 9 edges: the triangular prism and the triangular double pyramid; both achiral, so 2.
* 10 edges: there are two polyhedra, the pentagonal pyramid and a triangular prism with one square face
  cut along a diagonal. The second one is chiral (the diagonal can slope left or right), so the answer
  is 1.

      n:      6  7  8  9  10  11  12  13  14 ...       33
      a(n):   1  0  1  2   1   4   9  12  27 ... 17414512

(The entry starts at n = 1; the terms for fewer than 6 edges are 0.)

## 4. Why it's hard

Finding a shape's mirror symmetries is easy once you have the shape; the problem is having all the
shapes. The known values stopped at 33 edges. Counting them without building them needs a way to count
mirror-symmetric shapes directly, which is the new part of our method.

## 5. What we found

Forty-six new terms, from 34 to 79 edges:

**a(34) = 36,559,550**, **a(40) = 3,144,706,008**, …,
**a(79) = 27,080,101,035,925,607,425,280** (about 2.7·10^22).

The method is the one on the [main page](../layman.md): cut a mirror-symmetric polyhedron along its
mirror, count the halves with labels saying how each half meets the mirror, and add the
"turn and reflect" symmetries counted separately. This sequence is that mirror count, added up over all
shapes with n edges.

A striking pattern: mirror symmetry becomes very rare. With 10 edges, 1 shape in 2 is achiral. With
79 edges, about 1 in 3.7·10^15 is.

## 6. How we know it's right

* All 33 previously known terms come out exactly.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; the number with
  a mirror symmetry matched our formulas in every case.
* Swapping corners and faces turns an achiral polyhedron into another achiral one, so the counts must be
  symmetric under that swap. They are, even though our method treats corners and faces differently.
* The relation 2U − S ties this sequence to the two other edge counts, each checked against its own
  known terms (and it exposed a typo in one of them, see [A005645](layman_A005645.md)).

## 7. Why it matters (honestly)

This is an extension of a standard table, not a breakthrough.
* The table goes from 33 to 79 edges.
* It is the sequence that most directly tests the new mirror-counting idea: it is nothing but that count.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Achiral:** the same as its own mirror image (it has a "sense-reversing" symmetry).
- **Chiral:** different from its mirror image, like a hand.
