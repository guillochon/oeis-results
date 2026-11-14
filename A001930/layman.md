# A001930 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

Take a handful of points and decide which points are "close to" which. Mathematicians call a consistent
way of doing this a **topology**. On a finite set there is a simple way to picture one: draw arrows
between the points, with the rule that if there is an arrow from a to b and one from b to c, there must
be one from a to c. Two points with arrows both ways are "glued together": no open set can tell them
apart.

So a finite topology is a ranking with ties allowed, and with some pairs left unranked. Think of a
company chart where several people can share a desk (a tie), some people report to others, and some
pairs have nothing to do with each other.

## 2. Where this shows up

- Finite topological spaces are the simplest models of shapes, and they capture surprisingly much: every
  finite simplicial complex has a finite space with the same holes and loops.
- Rankings with ties and incomparable items appear in scheduling, version histories and preference
  data.
- The count is one of the classic "how many structures are there?" questions, alongside graphs, groups
  and partial orders.

## 3. What A001930 counts

The number of different topologies on n points, where two topologies are the same if one becomes the
other by renaming the points.

For 2 points there are 3:
- the two points are unrelated;
- one is below the other;
- the two are glued together.

The sequence begins 1, 1, 3, 9, 33, 139, 718, 4535, 35979, …

## 4. Why it's hard

The numbers grow explosively, and listing the structures one by one stops being possible around 16
points (5,637,294,117,525,695 topologies). That value, found by Brinkmann and McKay in 2005, stood as
the last known term for twenty years.

A topology in which no two points are glued is called a **partial order**. Counting topologies needs
the count of partial orders first, and that one had also been stuck at 16 points.

## 5. What we found

Earlier this month we counted the partial orders on 17 points (378,180,373,015,320,248). That
calculation worked by sorting partial orders according to their symmetries, and it left behind a table:
for every kind of symmetry, how many partial orders have it.

That table turns out to be all one needs. A topology is a partial order in which each point has been
replaced by a group of glued points, and the table says exactly how many ways there are to do that
without counting anything twice. The result:

**There are 457,507,251,208,437,773 topologies on 17 points.**

The same table answers a family of related questions with no new computation. Among them:
- the number of transitive relations on 16 and 17 points (**292,885,432,498,738,892,585** and
  **49,464,136,530,826,608,408,726**);
- the number of "prime" partial orders (those that cannot be built by plugging smaller ones into each
  other) for 13 to 17 points.

In all, 30 entries of the OEIS get new terms.

## 6. How we know it's right

- **Known values come out right.** The same program reproduces every previously known term of all 30
  entries, including the 17 known terms of A001930. This rules out errors in the formulas and, since each
  entry weighs the table differently, errors in individual table values for up to 16 points.
- **Whole numbers.** The formulas divide by numbers as large as 17! (about 3.6 × 10^14). Every result
  is a whole number, as it must be. A wrong table entry would almost surely leave a fraction.
- **Labeled versions.** When the points carry name tags the answers are known independently (to 19
  points); the same formulas reproduce them.
- **One honest caveat.** The 17-point values rest on our own 17-point count of partial orders. It passed
  the same kinds of checks and is now in the OEIS, but nobody else has recomputed it yet.

## 7. Why it matters (honestly)

This is a modest result with a nice moral, not a breakthrough.
- It extends a well-known reference sequence that had not moved since 2005.
- It shows that the symmetry table of partial orders, not just their count, is the valuable object: one
  table, thirty sequences.
- Several small entries that had 6 to 8 terms since the 1970s now have 18.
- The next step (18 points) is a clean, hard frontier: it needs partial orders on 18 points.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Topology (finite):** a set of points with a transitive "is below or glued to" relation.
- **Partial order:** a topology with no two points glued.
- **Transitive relation:** arrows such that a → b and b → c force a → c; a point may or may not have an
  arrow to itself.
- **Symmetry:** a renaming of the points that leaves the structure unchanged.
