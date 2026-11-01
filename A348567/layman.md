# A348567 in plain English

*How many sandpile patterns on a square grid look the same however you turn or flip the square? A formula that
answers it for every size, and a proof of a guess made last year.*

---

## 1. The sandpile

Take an n×n grid of cells and drop grains of sand on it. A cell that holds **four or more** grains is unstable:
it **topples**, sending one grain to each of its four neighbours. Cells on the edge lose the grains that would
go off the grid. One toppling can set off others, like an avalanche, until every cell holds 0 to 3 grains again.

This is the **abelian sandpile**, a toy model physicists introduced in 1987 for "self-organized criticality",
the way avalanches, earthquakes and forest fires come in all sizes.

## 2. Recurrent patterns

Keep dropping grains at random places and letting avalanches settle. Some stable patterns never come back once
you have left them; the ones that keep coming back are called **recurrent**. Remarkably, recurrent patterns
form a group: "add two patterns cell by cell, then let it topple" is a well-behaved addition. The group's
**identity** pattern, the one that changes nothing when you add it, is a famous fractal picture full of
nested squares and triangles.

## 3. What A348567 counts

a(n) is the number of recurrent patterns on the n×n grid that look exactly the same after any of the eight
symmetries of the square: the four rotations and the four reflections.

For n = 4 there are 12 of them; for n = 10 there are 2,876,576. Andrei Zabolotskii found the first eleven
values in 2021 and, in 2025, guessed a formula linking them to a neighbouring sequence (patterns that only
need to be symmetric left-right and top-bottom).

## 4. What is new

- **A formula for every n.** The count is a product of simple cosine expressions, one for each pair of odd
  numbers h ≤ k up to n:
  a(n) = product of (4 − 2cos(hπ/(n+1)) − 2cos(kπ/(n+1))).
- **The 2025 guess is now proved**, for every n.
- **New values**: a(11) through a(100). The last one has 635 digits.
- **A link to dominoes**: for even sizes, a(2m) is 2^m times A065072(m), a square root of the number of ways to
  tile a 2m×2m board with dominoes.

## 5. Why it works

A 2015 paper by Florescu, Morar, Perkinson, Salter and Xu showed that symmetric recurrent patterns are counted
by a determinant: that of the grid's "toppling matrix" restricted to symmetric patterns. A determinant is
the product of the matrix's eigenvalues. For the square grid the eigenvectors are waves, sine curves across and
down, and only the waves with odd frequencies h and k are symmetric. Making them symmetric across the diagonal
too pairs (h, k) with (k, h), so each pair counts once. Multiplying their eigenvalues gives the formula.

Compare it with the left-right and top-bottom symmetric count: there every ordered pair (h, k) counts. So that
count is a(n)² divided by the diagonal factors, and those multiply to an exact power of 2. That is the 2025 guess.

## 6. How we know the numbers are right

- For sizes up to 7, a computer tried every symmetric pattern and tested each one for recurrence. The counts
  match.
- For sizes up to 40, the determinant was computed exactly in whole numbers, with no rounding, and it equals
  the formula.
- For every even size up to 100, the values agree with the domino count A065072, which other people computed
  independently.
- All eleven values known before reproduce.

## 7. Credit

Laura Florescu's 2011 Reed College thesis had already linked these symmetric patterns, for even sizes, to
weighted domino tilings of a family of graphs studied by Lior Pachter. It gives no formula or values. The
determinant theorem is from Florescu, Morar, Perkinson, Salter and Xu (2015).
