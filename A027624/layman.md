# A027624 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

Picture a party where some pairs of guests have fallen out and refuse to stand next to each other.
You want to invite a group in which no two guests are feuding. How many such groups are there?
The empty group counts, and so does any single guest.

Mathematicians draw the guests as dots and the feuds as lines between dots. The drawing is a
**graph**, and a group with no line inside it is an **independent set**. Counting independent
sets is one of the most basic counting problems about graphs.

The graph here is the **hypercube**. Its dots are all the strings of n zeros and ones, and two
strings are joined by a line when they differ in exactly one position.
- For n = 2 it is a square: 00, 01, 11, 10 around the edges.
- For n = 3 it is an ordinary cube.
- For larger n it is the same pattern in higher dimensions.

## 2. Where this shows up

- **Statistical physics.** Independent sets are the "hard-core gas": particles that can't sit next
  to each other. Counting them is computing that model's partition function. The hypercube is a
  standard test case.
- **Combinatorics.** The growth rate on the hypercube was found by Korshunov and Sapozhenko in the
  1980s. Sapozhenko's proof uses the *graph container method*, now a standard tool that is still
  being refined (e.g. Jenssen and Perkins, 2019).
- **Coding theory and Boolean functions.** Strings of zeros and ones are codewords and inputs of
  logic circuits. An independent set in the hypercube is a set of codewords no two of which differ
  in a single bit.

## 3. What A027624 counts

**A027624(n)** is the number of independent sets in the n-dimensional hypercube.

Tiny example, n = 2 (the square 00 – 01 – 11 – 10 – 00):
- the empty set;
- the 4 single corners;
- the 2 diagonals, {00, 11} and {01, 10}.

That makes **7**. No set of three corners works, because some two of them would be neighbours.

The sequence starts 2, 3, 7, 35, 743, 254475, 19768832143. That last value, for n = 6, was the
end of the OEIS entry.

## 4. Why it looked hard

The numbers grow ferociously: each term is roughly the *square* of the one before. The programs in
the entry build the sets one corner at a time ("take this corner or not"), so their running time
is about as large as the answer. That was already tens of billions of steps for n = 6. The next
answer is about 4 billion times bigger, far beyond what that approach can reach.

## 5. What we found

> **A027624(7) = 78,685,477,899,897,082,403** (about 7.9 × 10¹⁹).
>
> **A027624(8) = 1,268,098,993,536,094,508,894,717,661,843,009,268,823** (about 1.3 × 10³⁹).

The key idea is to **stop listing sets and cut the cube into pieces**. A 7-dimensional cube can
be viewed as four 5-dimensional cubes arranged in a ring, each joined to its two neighbours.
An independent set of the big cube is then a choice of one independent set in each of the four
small cubes, where neighbouring choices never use the same position.

If you fix the choices in two opposite small cubes, the other two can be chosen separately.
Each of them just has to avoid the positions already used. So the count becomes a sum over pairs
of sets in a 5-dimensional cube, with no need to build any 7-dimensional set. We also use the
cube's symmetries: rotating or reflecting the small cube doesn't change the count, so only 288
essentially different first choices need checking. The whole computation takes under a second
on a home PC.

**One dimension up.** For n = 8 the same trick needs pairs of independent sets in a
6-dimensional cube, and there are about 20 billion of those sets, so about 10¹⁶ pairs even after
using the symmetries. That is too many to visit one by one. Two more ideas made it work.

- *Count the pairs instead of visiting them.* Many different pairs cover the same positions,
  and only the covered positions matter for the rest of the count. So we ask, for each set of
  positions Y, how many pairs cover exactly Y. The answer has a simple formula: look at the
  covered positions as a shape inside the cube. Each lone position can belong to the first set,
  the second, or both (3 ways). Each larger connected piece splits in exactly 2 ways. Multiply
  these up and you have the number of pairs. The sum now runs over shapes rather than pairs.
- *Do whole batches at once.* The shapes come in natural batches of 65,536 that differ only in
  one 16-corner slice of the cube. For such a batch the remaining count is a standard
  "sum over all subsets" transform, which computes all 65,536 values in one sweep. The cube's
  symmetries cut the number of batches to about 80 billion, organised around 1,228,158
  essentially different half-shapes.

The total work is still about 5 million billion small steps. A graphics card is built for exactly that
kind of repetitive arithmetic, so we wrote the inner loop for an ordinary gaming GPU (an RTX 3080).
Tuning it with the card's profiler brought the run down to **about 24 hours**. The answer is
bigger than a 128-bit integer, so the code adds in 192- and 256-bit pieces.

## 6. How we know it's right

- **Known answers reproduced.** The same program, applied to smaller cubes, gives back every
  earlier term (7, 35, 743, 254475, 19768832143). This rules out a wrong idea or a bug in the
  basic machinery.
- **Two ways of cutting.** Each earlier term also comes out of a second decomposition (a cube as
  two copies of a smaller cube, instead of four), which checks the helper routine separately.
- **Three runs, three levels of shortcut.** We computed the new value three times: using all the
  symmetries of the 5-dimensional cube (under a second), only the "flip a bit" symmetries
  (18 seconds), and no symmetry at all, checking all 65 billion pairs (10 minutes). All three
  agree to the last digit. This rules out a mistake in the symmetry shortcut.
- **Spot checks.** The fast helper that counts independent sets inside a region was compared with
  slow direct counting on 8000 random regions, with no disagreements.
- **The size is right.** Mathematicians know roughly how big these numbers must be. Our values
  are about 1.29 (n = 7) and 1.13 (n = 8) times the leading estimate, continuing the trend of the
  earlier terms (1.18, then 1.40) and closing in on 1 as the theory predicts. After the known
  first correction term is taken into account, what is left over is about 1.05 at n = 8.

For n = 8 specifically:
- **The shape-counting formula is checked on smaller cubes.** The same code, run on 3-, 4- and
  5-dimensional cubes instead of the 6-dimensional one, returns the known values for n = 5, 6
  and 7 exactly.
- **The graphics card is checked against the processor.** A separate, slower version of the
  count runs on the ordinary processor. The card's answers for 12 chosen shapes (including the
  empty shape and the full cube) and for 24 shapes picked at random after the run agree with it
  to the last digit, with and without the symmetry shortcut.
- **Nothing was skipped or double-counted.** The 1,228,158 essentially different shapes were
  listed by a separate program, their multiplicities add up to exactly 4,294,967,296 = 2³², and
  every one of them appears exactly once in the result file.
- **What has not been done.** A second full 24-hour run without the symmetry shortcut would be a
  completely independent recomputation of the total. It is planned but has not been run yet.

## 7. Why it matters (honestly)

This is a small result, not a breakthrough.

- It extends a classic, well-studied sequence by two terms. The sequence had been stuck at n = 6,
  and the hypercube is a standard benchmark for the asymptotic theory.
- Exact values let researchers test how quickly the known formulas for the growth rate
  kick in. For example, one can now see how much the correction terms matter at n = 7 and 8.
- It is a nice example of old terms that *looked* hard only because of the method. Cutting
  the cube into pieces turned an impossible listing job into a one-second calculation, and
  counting pairs by their footprint turned the next term into a day on a gaming GPU.
- The next term, n = 9, is about 10⁷⁸ and remains out of reach. The same trick one level up would
  need every subset of a 7-dimensional cube's 128 corners, which is hopeless, so it needs a
  genuinely new idea.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Graph:** dots (vertices) joined by lines (edges).
- **Independent set:** a set of dots with no line between any two of them. The empty set counts.
- **Hypercube Qₙ:** the graph on all strings of n zeros and ones, with a line between strings that
  differ in one position.
- **Symmetry:** a rearrangement of the cube (flipping a bit everywhere, or swapping two positions)
  that maps lines to lines.
