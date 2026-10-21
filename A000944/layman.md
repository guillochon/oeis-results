# A000944 in plain English

*How many different shapes of polyhedra have 19 or more corners?*

---

## 1. The idea

A polyhedron is a solid with flat faces: a cube, a pyramid, a soccer ball. Two polyhedra have the same
*shape* here if they are built the same way, whatever the lengths and angles. Same number of corners,
connected by edges in the same pattern, and the same faces. Squash a cube into a slanted box and it is
still "the cube". Cut one corner off and you get a new shape.

So the question is: how many genuinely different ways are there to build a polyhedron with n corners?

## 2. Where this shows up

* Chemistry and crystallography: the shapes of molecules and crystal cages are polyhedra.
* Graph theory: by Steinitz's theorem (1922) these shapes are exactly the networks that can be drawn in
  the plane without crossings and stay connected after removing any two points. Counting polyhedra is
  counting those networks.
* Computer graphics and geometry use catalogues of small polyhedra as test cases.

## 3. What A000944 counts

The number of polyhedron shapes with n corners. With 4 corners there is only the triangular pyramid
(the tetrahedron). With 5 corners there are two: the square pyramid and the triangular double pyramid.
With 6 corners there are 7, among them the octahedron and the triangular prism.

    n:      4  5  6   7    8     9     10      11 ...      18
    a(n):   1  2  7  34  257  2606  32300  440564 ... 107854282197058

## 4. Why it's hard

The numbers grow by a factor of about 17 for each extra corner. The known values were found by
computer programs that *build every shape once* and count them. The best of these (plantri, by Gunnar
Brinkmann and Brendan McKay) reached 18 corners in 2006: about 10^14 shapes. For 19 corners there are
about 2·10^15, too many to build one by one, and the sequence had stopped there.

## 5. What we found

**a(19) = 1,837,856,160,118,250**, **a(20) = 31,696,106,245,572,086**,
**a(21) = 552,533,560,690,636,208**, and so on up to 40 corners:
a(40) = 680,349,069,787,813,513,725,946,584,317,522,497,406,424 (about 7·10^41).

Instead of building the shapes, we counted them with formulas. Counting shapes "with a labelled
starting edge" is easy: there is a 1960s formula for that. The difficulty is the symmetric shapes.
A shape with symmetries looks the same from several starting edges, so it would be counted several
times. A classical trick (Burnside's lemma) fixes this, if you can also count how many shapes have
each kind of symmetry: a rotation, a mirror, or a "turn and reflect". We found ways to count each
kind. For mirrors, the new idea is to cut a symmetric polyhedron along its mirror and count the halves,
together with labels saying how each half meets the mirror. For "turn and reflect" symmetries, folding
the polyhedron onto itself gives a map on the projective plane (a sphere with opposite points glued),
and those maps can be counted with a 1960s-style equation.

## 6. How we know it's right

* For up to 11 corners we built every polyhedron by computer and listed all its symmetries. Every
  symmetry count from the formulas matched, number by number. This rules out a wrong formula for any
  kind of symmetry.
* The same formulas give all the known values, for 4 to 18 corners, exactly. They also match
  hundreds of published finer counts: by corners and faces together, by edges, and the count of
  mirror-symmetric polyhedra by edges.
* Swapping corners and faces turns every polyhedron into another one (cube ↔ octahedron), so the
  counts must be symmetric under that swap. Our mirror method treats corners and faces differently,
  and the counts still came out symmetric.
* The answer is a sum of fractions, and it came out a whole number in every case, as it must.
* Mathematicians proved decades ago how fast these numbers must grow in the long run. Our new numbers
  approach that growth rate smoothly, exactly as the theory predicts.

## 7. Why it matters (honestly)

This is one more term of a classic sequence. It is not a breakthrough, but this sequence is a central
one.
* A000944 is one of the OEIS's "core" sequences; these are its first new terms in 20 years, and they
  more than double its length: from 18 to 40 corners.
* The method needs no enumeration, so it is not limited by how many shapes exist. All 22 new terms took
  under two hours on a desktop PC.
* The mirror-counting idea is new as far as we know and should apply to other kinds of maps.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Polyhedron shape (combinatorial type):** what is left of a polyhedron when you forget lengths and
  angles and keep only which corners, edges and faces touch.
- **Burnside's lemma:** a rule for counting things "up to symmetry" by averaging over the symmetries.
