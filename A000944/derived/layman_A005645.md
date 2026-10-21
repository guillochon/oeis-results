# A005645 in plain English

*How many polyhedra have n edges, if a shape and its mirror image count as two?*

---

## 1. The idea

Your left and right hands are mirror images, but no turning will make one fit into the other's glove.
Some polyhedra are like that too: they are *chiral*. Others, like the cube, look the same in a mirror.

As on the [main page](../layman.md), a polyhedron's shape is the pattern of its corners, edges and faces,
forgetting lengths and angles. This sequence counts shapes by number of edges, with one change: a chiral
shape and its mirror image count as *two* different shapes. You may turn a polyhedron around, but not
reflect it. Mathematicians call this counting "sensed" maps.

## 2. Where this shows up

* Chemistry: mirror-image molecules (left- and right-handed versions) can behave very differently, so
  telling them apart matters.
* Counting theory: the "mirror images separate" count is the one that classical map-counting formulas
  handle most naturally; the count with mirror images identified is built from it.

## 3. What A005645 counts

The number of polyhedra with n edges, mirror images counted separately.

The smallest chiral polyhedron has 10 edges. Take a triangular prism and cut one of its square faces
into two triangles along a diagonal. The diagonal can slope "left" or "right", and no turning of the
solid changes one into the other. So with 10 edges there are 3 shapes in this count: the pentagonal
pyramid (its own mirror image) and the twisted prism in both versions. The ordinary count of
[A002840](layman_A002840.md) has 2.

    n:      6  7  8  9  10  11  12  13  14 ...                50
    a(n):   1  0  1  2   3   4  15  32  89 ... 3576732121617644389508

## 4. Why it's hard

The numbers grow about 3.6-fold per extra edge. The entry's terms up to 50 edges came from earlier
calculations; beyond that nobody had gone.

## 5. What we found

**A correction.** The entry says a(34) = 3,323,25**8**,053,528. The right value is
**a(34) = 3,323,268,053,528**: one digit is wrong.

**Twenty-nine new terms**, from 51 to 79 edges: **a(51) = 13,330,027,831,373,735,087,282**, …,
**a(79) = 202,513,714,771,293,355,444,916,723,643,665,655,888** (about 2·10^38).

They come from the [A000944](../layman.md) computation, which counts polyhedra by corners and faces using
formulas and Burnside's lemma, separately for "turning only" and "turning or reflecting". This sequence
adds up the "turning only" counts over all shapes with n edges.

## 6. How we know it's right

* Of the 45 previously known terms (6 to 50 edges), 44 come out exactly; the other is a(34).
* **Why a(34) is a typo, not our mistake.** Every shape is either its own mirror image or one of a mirror
  pair. A mirror pair counts twice here and once in the ordinary count; a mirror-symmetric shape counts
  once in both. So

      (this count) + (mirror-symmetric shapes) = 2 × (ordinary count).

  Check with 10 edges: 3 + 1 = 2 × 2. For 34 edges the ordinary count A002840(34) = 1,661,652,306,539 is
  known independently, from a computer program that built every polyhedron. Our mirror-symmetric count
  is 36,559,550, and 3,323,268,053,528 + 36,559,550 = 3,323,304,613,078 = 2 × 1,661,652,306,539 exactly.
  With the stored value the left side is 10,000,000 short. The discrepancy is exactly one digit in one
  place, the signature of a copying error.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; every count from
  the formulas matched.
* Swapping corners and faces must leave the counts unchanged, and it does.

## 7. Why it matters (honestly)

This is a small result: a standard table extended, and one typo fixed.
* The correction matters most: a wrong digit in a reference table can quietly break anyone's check that
  relies on it. It had sat in both the entry's terms and its b-file.
* The table now runs from 50 to 79 edges.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Chiral:** different from its mirror image, like a hand.
- **Sensed map:** a shape counted up to turning only, so mirror images are distinct.
- **b-file:** the long list of terms attached to an OEIS entry.
