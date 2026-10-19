# A005162 in plain English

*Counting a very symmetric kind of "alternating sign" grid: every size up to 62 (the OEIS stopped at 19), and the
first new even-size counts since 1992.*

---

## 1. The idea
Fill a square grid with the numbers 1, 0 and −1 so that:
- every row and every column adds up to 1, and
- along any row or column, if you skip the zeros, the nonzero numbers alternate +1, −1, +1, …, starting and
  ending with +1.

These grids are called **alternating sign matrices** (ASMs). If you never use −1, they are permutation matrices
(one 1 in each row and column), so ASMs are a generalisation of permutations. In the 1980s it was found, to
everyone's surprise, that the number of n×n ASMs is given by a beautiful product formula
(1, 2, 7, 42, 429, …). Proving it took over a decade and became a famous story in combinatorics.

## 2. Where this shows up
- **Statistical physics:** ASMs correspond exactly to configurations of "square ice", a lattice model of water
  molecules (the six-vertex model). The famous proof of the counting formula came from this physics.
- **Algebra:** ASMs appear in the study of determinants (Robbins and Rumsey's λ-determinant) and in
  representation theory.
- **Plane partitions:** ASMs are counted by the same numbers as some families of 3-D stacks of cubes, and nobody
  yet has a simple direct explanation of why.

## 3. What A005162 counts
Mathematicians then asked how many ASMs have each kind of symmetry of the square: mirror symmetry, rotation, and
so on. A005162 counts the ASMs that are **symmetric about both diagonals**: reflecting the grid in the main
diagonal, or in the other diagonal, leaves it unchanged.

Example, n = 3. The 3 such grids are

```
1 0 0     0 0 1     0 1 0
0 1 0     0 1 0     1 -1 1
0 0 1     1 0 0     0 1 0
```

The counts start 1, 2, 3, 8, 15, 52, 126, 568, 1782, …

## 4. Why it's hard
For **odd** sizes there is again a neat product formula. It was guessed in the 1980s and finally proved in 2017.
For **even** sizes nobody has found a formula, and experts believe there is none: the numbers have huge, random-
looking prime factors. So the only way to get even-size values is to count them. The numbers grow like
1.3^(n²/4), and a(62) has 117 digits, so listing the grids one by one is hopeless. Until now the published values
stopped at size 24 (a 1992 table), and the OEIS entry stopped at size 19.

## 5. What we found
We computed **every value up to size 62**. The 1992 table already had sizes 20, 21, 22 and 24, and our values
agree with it. The values for sizes 23 and 25 to 62 are new, including the first new even sizes, 26, 28, …, 62.
For example,

- **a(26) = 19990299994021601231552**
- **a(62)**, a 117-digit number.

The key idea: because of the two symmetries, the whole grid is fixed by one triangular quarter of it. The rules
"rows and columns add up correctly" become a square-ice model on that triangle, with simple rules where the
triangle meets the diagonals. Then we sweep across the triangle column by column. At each moment we only need to
remember about n/2 yes/no facts (whether each half-finished row has a running sum of 0 or 1), so the computer
tracks about 2^(n/2) cases, not the astronomically many grids. Size 62 took 7 minutes on a desktop PC.

## 6. How we know it's right
- **It reproduces all 19 values already in the OEIS.**
- **Odd sizes check against the proved formula.** The same program, with a one-cell difference at the centre,
  gives the odd sizes, and all 31 odd values up to 61 match the 2017 formula exactly. A bug in the method would
  almost certainly break these.
- **It matches the 1992 table** at sizes 20, 21, 22 and 24 (computed independently).
- **Three separate programs** (plain Python, numpy, and a fast Rust version with two different number formats)
  agree wherever they overlap. Each big number was computed modulo many primes and reassembled, and one extra
  prime was kept back to confirm the reassembly.
- **The growth is right:** the values approach the proved growth rate 3√3/4.

Along the way we noticed that the 1992 table's size-9 value (1783) is a misprint for 1782.

## 7. Why it matters (honestly)
This is a solid data result, not a breakthrough.
- It extends an OEIS entry that had been stuck at size 19 since the 1990s, to size 62.
- The even sizes are exactly the open case. Anyone looking for a formula, an asymptotic expansion, or patterns in
  the prime factors now has 39 more values than were published before (43 more than the OEIS had).
- The method (a column sweep with a stack of open rows) is short and reusable for the related symmetry classes.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number sequences.
- **Alternating sign matrix (ASM):** a grid of 1, 0, −1 with row and column sums 1 and alternating nonzero signs.
- **Square ice / six-vertex model:** a physics model of ice on a grid; its configurations correspond to ASMs.
- **Product formula:** a formula that writes a count as a product of simple factors, such as factorials.
- **Chinese remainder theorem:** a way to rebuild a huge number from its remainders modulo several primes.
