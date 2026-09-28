# A114601 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

Take a square grid of numbers, n rows by n columns, that is **symmetric**: the entry in row i,
column j equals the entry in row j, column i. Put a 2 in every spot on the main diagonal, and
fill every other spot with −1, 0 or 1.

Some of these grids are **positive definite**. That is the matrix version of "positive". Such a
matrix describes a genuine geometry: there are n arrows (vectors) in space whose lengths and
angles are recorded by the matrix. Here, the number in row i, column j is the "dot product" of
arrows i and j. The 2s on the diagonal say every arrow has the same length (√2). The ±1s and 0s
say two arrows meet at 60°, 120° or 90°.

**A114601(n)** counts the positive definite grids of size n.

Tiny example, n = 2: the grid is `[[2, x], [x, 2]]` with x ∈ {−1, 0, 1}. All three choices are
positive definite (two arrows at 120°, 90° or 60°), so a(2) = **3**.

## 2. Where this shows up

- **Crystals and sphere packings.** A set of arrows like this generates a *lattice*, a regular
  grid of points in space, like atoms in a crystal. The lattices you get from arrows of length √2
  are exactly the famous **root lattices**. There are only five families of them, called A, D, E6,
  E7 and E8. E8 gives the densest way to pack spheres in eight dimensions.
- **Lie theory and physics.** The same five families classify the "simply-laced" symmetry groups
  that appear throughout mathematics and theoretical physics. The E8 lattice turns up in string
  theory.
- **Matrix counting.** The entry is a sibling of A086215 (all positive definite matrices with
  entries −1, 0, 1) and A085657 (the same with only 0s and 1s off the diagonal).

## 3. What was known

The entry listed seven terms: 1, 3, 23, 393, 13089, 737595, 58969079. They came from a program
that tries every grid and tests each one. There are 3^(n(n−1)/2) grids to try, which already
reaches about 10^17 for n = 8. The entry was marked *hard*, and nobody had a formula.

## 4. What we found

**A formula that gives every term.** The key fact is the lattice connection from section 2. Every
grid we count is the "angle table" of a set of arrows that build one of the root lattices A, D, E.
Classifying the lattices then classifies the grids:

- A grid splits into independent blocks (arrows at right angles to everything outside their
  block). Each block builds a single A, D or E lattice. Counting whole grids from blocks is a
  standard trick (the *exponential formula*).
- For the **A** family, the ways of choosing the arrows correspond to *trees*: networks with no
  loops. Trees on labelled points were counted by Cayley in 1889.
- For the **D** family, they correspond to networks with **exactly one loop**, where the loop has
  an odd number of "+" links. Those are also classical to count.
- The **E** family has just three members, E6, E7 and E8. Each gives one fixed number, which we
  computed by computer, two different ways.

Put together, this gives a single formula. It involves the "tree function" of combinatorics plus
three constants: 355104, 43733504 and 4053973888. From it, a computer can produce hundreds of
terms in seconds:

> a(8) = **4,925,177,553**, a(9) = **63,545,326,593**, a(10) = 1,772,133,404,723, …

The formula also shows how fast the numbers grow: roughly n! × (2e)ⁿ, with an exact constant in
front.

## 5. How we know it's right

- **Proof.** The formula is proved: a written proof, built on a classical theorem of Witt (1941)
  that lists all root lattices. The two graph descriptions (trees for A, one-loop networks for D)
  are also checked by the **Lean** proof assistant, which verifies every logical step by computer.
- **Old terms reproduced.** The formula gives back all seven known terms.
- **Brute force, twice.** We wrote a fast program that really does build every positive definite
  grid, one row at a time, with exact whole-number arithmetic. It confirms a(8) and a(9). For n = 9
  it took about 6 hours on 12 processor cores. It also confirms the finer prediction of *which*
  lattices appear: at size 9 only the A9 and D9 types can form a single block, and their counts
  (2,560,000,000 and 13,545,271,296) match the formula exactly.
- **The E numbers, two ways.** Each E constant was computed once by counting grids and once by
  working inside the actual E6, E7 and E8 lattices. The second count was checked against a
  classical identity (Cauchy–Binet) that must hold exactly.

## 6. A curious bump

The terms don't grow perfectly smoothly. Around n = 16 the sequence takes an unusually big jump,
then a smaller one. That is the E8 lattice showing off: at n = 16 you can build two independent
copies of E8, and that alone accounts for about 40% of a(16). The effect fades for larger n.

## 7. Why it matters (honestly)

This is a small piece of mathematics, not a breakthrough. It still has real value:

- **A "hard" sequence becomes an easy one.** The entry was marked *hard*: its terms could only be
  found by trying an astronomically growing number of grids. Now there is an exact, proved formula,
  and any term can be computed in moments. That is the best outcome an OEIS entry can hope for.
- **Numbers are now possible for sizes that were out of reach.** Brute force stalled at n = 7, and
  our own brute-force check needed 6 hours for n = 9. The formula gives every term up to n = 200,
  and more if needed, plus exactly how fast the numbers grow.
- **A classical theorem at work.** Witt's classification of root lattices, from the 1940s, does
  almost all of the counting. It is a nice example of deep structure turning an ugly search into
  a clean formula. It also explains the odd bump at n = 16, which would be baffling without E8.
- **A method, not just an answer.** The same idea should apply to the sibling entries A086215
  (all positive definite −1/0/1 matrices) and A085657 (only 0s and 1s off the diagonal). We plan
  to try A086215 next.
- **Honest limits.** Every ingredient is classical, so this is not research-level new mathematics.
  The new part is putting the ingredients together for this count, and the new terms. The three E
  constants rest on computer counts, checked two independent ways, rather than on a proof by hand.

## 8. Where to look

- The OEIS entries: [A114601](https://oeis.org/A114601), and its siblings
  [A086215](https://oeis.org/A086215) and [A085657](https://oeis.org/A085657).
- [A320064](https://oeis.org/A320064): the "reflectable bases" of type D, the same count as our
  D-family blocks. It comes from S. Azam, M. B. Soltani, M. Tomie and Y. Yoshii, *A graph
  theoretical classification for reflectable bases*, 2019.
- E. Witt, *Spiegelungsgruppen und Aufzählung halbeinfacher Liescher Ringe*, 1941: the
  classification of root lattices.
- J. H. Conway and N. J. A. Sloane, *Sphere Packings, Lattices and Groups*, Chapter 4: a friendly
  account of the root lattices A, D, E.
- Our full proof, with the Lean formalization and the computer checks: [`proofs.md`](proofs.md).

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Symmetric matrix:** a square grid of numbers that is unchanged by flipping it across its main
  diagonal.
- **Positive definite:** the matrix is the table of dot products of some set of independent arrows.
- **Lattice:** all the points you can reach by adding and subtracting a fixed set of arrows.
- **Root lattice:** a lattice built from arrows of length √2 (the "roots"). There are only the
  types A, D, E6, E7, E8 and their combinations.
- **Tree:** a network of points and links with no loops.
