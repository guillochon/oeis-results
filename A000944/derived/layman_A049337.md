# A049337 in plain English

*A table of polyhedra by corners and edges, and why every entry in it was in the wrong place.*

---

## 1. The idea

As on the [main page](../layman.md), a polyhedron's shape is the pattern of its corners, edges and
faces, forgetting lengths and angles. The main page counts shapes by corners. This entry is a finer
table: row n lists, for each possible number of edges k, how many shapes have n corners and k edges.

A network with n corners can have at most n(n − 1)/2 links (every pair joined), so row n has an entry
for k = 0, 1, 2, …, n(n − 1)/2. Most are 0: a polyhedron with n corners always has between 3n/2 and
3n − 6 edges.

## 2. Where this shows up

* Graph theory: by Steinitz's theorem these shapes are the 3-connected planar graphs, so this is the
  standard "by nodes and edges" table for them.
* Each row adds up to the main count A000944(n), and the last nonzero entry of row n counts triangle-only
  polyhedra ([A000109](layman_A000109.md)).

## 3. What A049337 counts

The rows, starting from n = 1 corner. Rows 1 to 3 have 1, 2 and 4 entries, all 0 (no polyhedra).

* Row 4 (k = 0..6): 0, 0, 0, 0, 0, 0, **1**. The tetrahedron has 4 corners and 6 edges.
* Row 5 (k = 0..10): only k = 8 (square pyramid) and k = 9 (triangular double pyramid) are 1.
* Row 6: 1 shape with 9 edges (the triangular prism), then 2, 2 and 2 with 10, 11 and 12 edges:
  7 in all, as on the main page.

## 4. What was wrong

The stored table had lost its first three entries, all zeros. Rows 1 to 3 should hold 1 + 2 + 4 = 7
entries, but held 4. Since the table is read as one long list, every later number sat three places
early. Row 4 read 0, 0, 0, 1, 0, 0, 0, which puts the tetrahedron at 3 edges, and every row after that
was shifted the same way. (A second table in the entry, written the other way round, was right.)

## 5. What we found

**The correction:** put three zeros back at the start. Then all 469 stored terms (rows 1 to 14) agree
with our computation, number for number.

**The extension:** rows 15 to 39, for 9,919 entries in all. For example, row 15 runs from 1,404 shapes
with 23 edges up to 2,406,841 with 39 edges; row 39 runs from 75,441,790,558,926 shapes with 59 edges
to **16,583,511,070,523,241,238,223,171,006** with 111 edges (the triangle-only ones), peaking at
**4,231,248,875,578,461,331,552,625,778,501,401,307,147** with 90 edges.

These are entries of the corners-and-faces table from the [A000944](../layman.md) computation
(formulas and Burnside's lemma, no shapes built), relabelled: by Euler's rule a shape with n corners and
k edges has k − n + 2 faces.

## 6. How we know it's right

* With the three zeros restored, all 469 previously stored terms match. A shift by exactly three
  places that makes every one of them agree is not a coincidence.
* Every row adds up to A000944(n), which matches the long-known values up to 18 corners, and the
  table's last entries in each row match
  [A000109](layman_A000109.md) and other published counts.
* For up to 11 corners we built every polyhedron by computer and listed its symmetries; every count from
  the formulas matched.
* Swapping corners and faces must leave the underlying table unchanged, and it does.

## 7. Why it matters (honestly)

This is a small result. The new rows are a by-product; the main value is the fix.
* Anyone reading the old table literally would have got wrong answers for every shape (a tetrahedron
  with 3 edges). Restoring three zeros makes the entry match its own definition.
* The table now runs from 14 to 39 corners.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Triangle read by rows:** a table stored as one long list, row after row, so a missing entry shifts
  everything after it.
- **Euler's rule:** for any polyhedron, corners − edges + faces = 2.
