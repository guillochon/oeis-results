# A005162: n×n alternating sign matrices symmetric about both diagonals (DASASMs), n ≤ 62

An **alternating sign matrix** (ASM) is a square matrix with entries −1, 0, 1 in which every row and every column
sums to 1 and the nonzero entries alternate in sign along each row and column. A **DASASM** is an ASM that is
symmetric about both diagonals: A = Aᵀ and A(i, j) = A(n+1−j, n+1−i). A005162(n) is the number of n×n DASASMs.

**Result.** a(n) for every n ≤ 62 (`b005162.txt`; a(62) has 117 digits). Before this, the OEIS had a(1)–a(19).
Bousquet-Mélou and Habsieger's table (1995, below) also had a(20)–a(22) and a(24), and our values agree with
them. a(23) and a(25)–a(62) are new. For odd n they also follow from the product formula proved by
Behrend, Fischer and Konvalinka. **The even terms a(26), a(28), …, a(62) are new: no formula is known or
conjectured for even n.**

```
a(20) = 55203533023856
a(21) = 828751754742975
a(22) = 23322499520009200
a(23) = 454739496669274500
a(24) = 16623477112752060736
a(25) = 420972227408592675000
a(26) = 19990299994021601231552
a(27) = 657522745057190417409000
a(28) = 40558893635955279908257088
a(29) = 1732789066323343611643088400
a(30) = 138844899597490901146947621952
...    (n ≤ 62 in b005162.txt)
```

## Background
- DASASMs are one of the symmetry classes of ASMs studied by Robbins (*Symmetry classes of alternating sign
  matrices*, arXiv:math/0008045). Odd order has the product formula
  |DASASM(2m+1)| = ∏_{i=0}^{m} (3i)!/(m+i)!. Robbins conjectured it, and Behrend, Fischer and Konvalinka proved it
  (*Diagonally and antidiagonally symmetric alternating sign matrices of odd order*, Adv. Math. 315 (2017)
  324–365, arXiv:1512.06030).
- **Even order** has no known or conjectured formula. Behrend, Fischer and Koutschan (*Diagonally symmetric
  alternating sign matrices*, arXiv:2309.08446) list it, with diagonally symmetric ASMs (A005163) and odd totally
  symmetric ASMs (A005164), as the three classes where prime factorizations seem to rule out a product formula.
  They prove the growth rate lim |DASASM(n)|^{4/n²} = 3√3/4.
- Data: M. Bousquet-Mélou and L. Habsieger, *Sur les matrices à signes alternants*, Discrete Math. 139 (1995)
  57–72, Table 1 (column A₄(n)): n ≤ 22 and n = 24. (Their n = 9 entry, 1 783, is a misprint for 1782: the OEIS,
  the product formula and our count all give 1782.) We found no published values beyond n = 24.

## Method: a six-vertex model on the fundamental triangle, swept column by column

**Reduction to a triangle.** A DASASM is determined by its entries on
T = {(i, j) : i ≤ j ≤ n+1−i}, the top quarter between the two diagonals. For i ≤ ⌈n/2⌉, row i of A reads

  S_i = (A(1,i), …, A(i−1,i) | A(i,i), …, A(i,n+1−i) | A(i−1,n+1−i), …, A(1,n+1−i)).

That is column i of T downward (since A(i, j) = A(j, i)), then row i of T, then column n+1−i of T upward (since
A(i, j) = A(n+1−j, n+1−i)). Every other line of A is one of these or a reversal of one: columns equal rows by
A = Aᵀ, and row i > n/2 is row n+1−i read backwards. Reversal keeps the ASM line condition, because suffix sums
are 1 minus prefix sums. So **A is a DASASM iff every S_i has partial sums in {0, 1} and total 1.**

**The model.** Each interior cell of T (off both diagonals) lies on exactly two of these sequences: row r of T,
and the column part of the hook through column c. Write h for the running row sum entering a cell from the left
and v for the running column sum (top-down) entering it from above. The partial-sum conditions become:

- **interior cell:** entry x with h+x, v+x ∈ {0, 1}. So (h, v) = (0, 0) → (0, 0) or (1, 1); (1, 1) → (1, 1) or
  (0, 0); (0, 1) and (1, 0) pass through unchanged (x = 0). This is the usual square-ice rule.
- **diagonal cell (j, j)**, on the main diagonal: the column-j sum v so far is the prefix of S_j. Row j starts
  here with running sum v + A(j, j), which may be 0 or 1 whatever v is (the entry is then determined and lies in
  {−1, 0, 1}). So **the new row sum is a free bit** and v is forgotten.
- **antidiagonal cell (n+1−j, j):** row n+1−j ends here, and so does the column-j sum v of the cells above. The
  total h + x + v = 1 **forces** x = 1 − h − v, and it is always valid: x ∈ {−1, 0, 1} and the running sum 1 − v is
  in {0, 1}. The upward part of S_{n+1−j} needs top-down partial sums of column j in {0, 1}; the interior rule
  already ensures that.
- **centre cell** (n odd, j = (n+1)/2): S_j is v, then A(j, j), then v again, so x = 1 − 2v. Always valid.

So |DASASM(n)| is the number of ice configurations on T with free turns on the diagonal and forced ends on the
antidiagonal.

**The sweep.** Process the columns of T from left to right, each from the top down. Row j opens at the diagonal
cell of column j and closes at the antidiagonal cell of column n+1−j. The open rows therefore form a stack, and
the state between columns is just the running sums of the open rows: at most ⌊n/2⌋ bits, plus the one vertical
bit while a column is in progress. A row-by-row sweep would need about n bits. For even n the middle row n/2 opens
at column n/2 and closes at column n/2+1, with no cell in between touching it, so it is never stored (a factor 2).
Peak storage is 2^⌈n/2⌉ counters. Each cell update is in place. An interior cell adds the (0,0) and (1,1)
entries of each pair of states and writes the sum to both. Opening a row duplicates the array, and closing one
sums out two bits.

Counts are taken modulo primes and combined by the Chinese remainder theorem. Enough primes are used to exceed
2·|DASASM(2⌊n/2⌋+1)|, an upper bound because |DASASM(n)| ≤ |DASASM(n+1)| (Behrend–Fischer–Koutschan), and one
**spare prime** is used as a check.

**Worked example, n = 4.** T has cells (1,1), (1,2), (1,3), (1,4), (2,2), (2,3). Row 1 opens at (1,1) with a free
bit r₁. Cell (1,2) is interior with v = 0. If r₁ = 0 it gives two outcomes, h = 0 or h = 1 (the latter with
v = 1 going down). If r₁ = 1 it gives one outcome, h = 1. So there are 3 ways so far. Cell (2,2) opens row 2
with a free bit (the middle row: factor 2). Cell (1,3) is interior with v = 0: an incoming h = 0 gives 2
outcomes and h = 1 gives 1, for 2 + 1 + 1 = 4 ways. Cells (2,3) and (1,4) are antidiagonal ends (forced).
Total 4 · 2 = **8** = a(4).

## Checks
- **Reproduces all 19 terms** previously in the OEIS.
- **Odd n: all 31 odd values n ≤ 61 equal the proved product formula exactly.** The program has no special case
  for odd n except the single centre cell, so this checks the same transfer matrix that gives the even terms.
- **Agrees with Bousquet-Mélou–Habsieger** at n = 20, 21, 22, 24 (an independent computation from 1992).
- **Three implementations agree.** `dasasm.py` (exact Python integers, a dict of states) agrees with
  `dasasm_fast.py` (numpy, 62-bit primes) for n ≤ 30. `dasasm_fast.py` agrees with the Rust program for n ≤ 50,
  in both its builds (u32 residues mod 31-bit primes, and u16 residues mod 16-bit primes).
- The spare prime agrees for every n.
- **Monotone:** a(n) ≤ a(n+1) for all n < 62, as it must be (Behrend–Fischer–Koutschan).
- **Growth:** a(n)^{4/n²} = 1.372 (n = 20), 1.347 (30), 1.335 (40), 1.328 (50), 1.322 (62), decreasing slowly
  towards the proved limit 3√3/4 = 1.299.
- **No product formula, again:** the even terms have huge prime factors. a(24) = 2⁶·19·149·163·562878191033,
  a(26) = 2⁶·19·53·310177197027395749, a(28) = 2⁶·7·11·193·42644015413619625097. This agrees with
  Behrend–Fischer–Koutschan's observation.

## Programs and reproduction
- `dasasm.py NMAX`: exact counts (dict of frontier states), n ≤ ~32. About 60 lines, and the clearest statement
  of the model. Also has the odd-order product formula.
- `dasasm_fast.py NMIN NMAX`: numpy, mod primes + CRT. n = 50 takes 98 s.
- `rust/` + `run_rust.py NMIN NMAX [--w16]`: the in-place multithreaded version.
  `cargo build --release` gives u32 storage (primes < 2^31). `CARGO_TARGET_DIR=target16 cargo build --release
  --features w16` gives u16 storage (primes < 2^16, half the memory). `run_rust.py` does the CRT, the spare-prime
  check and the comparisons above. Timings with 10 threads on an i5-12600K: n = 50 in 4 s; n = 56 in 36 s (u32,
  12 primes); n = 62 in 407 s (u16, 28 primes, 4 GB). n = 63–64 would need 8 GB in u16.
- `results/`: `a005162_1_62.txt` (all values), `exact_1_30.txt`, `fast_1_50.txt`, `rust.txt` (one line per run:
  n, a(n), primes, seconds; n = 40–44 appear twice, once from each build), `rust16_57_62.log`. The u32 run for
  n = 53–56 is in `rust.txt`; it was stopped after n = 56 to switch to the u16 build for memory.

## Working log
- 2026-10-19: chosen as the highest-significance open target. Literature: Valency and OpenAlex searches for
  DASASM / "antidiagonally symmetric", and the TeX sources of arXiv:1512.06030, 1611.03823 and 2309.08446.
  The only data beyond the OEIS is Bousquet-Mélou–Habsieger's table, n ≤ 24.
- Column-sweep model derived. `dasasm.py` reproduced a(1..19) and the odd formula to n = 21 on its first run.
  numpy to n = 50, then Rust (in place, unstored middle row, u16 build) to n = 62.
