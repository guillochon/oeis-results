# A157615 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea: a walk that must turn at every step

Take an n × n chessboard and walk from square to square. Each step goes to a neighbouring square
(left, right, up or down), and you may never visit a square twice. There's one twist: **the steps
must alternate** between horizontal and vertical. Right, then up or down, then left or right, then
up or down, and so on. You can never take two horizontal steps in a row, so the walk turns 90° at
every single square.

How many squares can such a walk cover? That maximum is **A157615(n)**.

Here's a best possible walk on a 4 × 4 board. S is the start, E is the end, `#` marks the other
visited squares, and `.` the two squares it misses. It covers 14 of the 16 squares:

    S # # E
    # # # #
    # # # #
    . # # .

## 2. Where this shows up

- **Puzzles:** it's a close cousin of the knight's tour and of "snake" puzzles. It came from an
  idea of Leroy Quet discussed on the SeqFan mailing list.
- **Path problems in grids:** questions like "how long can a restricted self-avoiding walk be?"
  come up in computer science (routing on chips, robot motion with limited turning) and in the
  physics of polymers, where self-avoiding walks model long molecules.
- **A test case for proof methods:** it's simple to state but surprisingly stubborn, which makes
  it a good small example of how to prove that something is *impossible*.

## 3. What A157615 counts

For boards of size 1, 2, 3, … the maximum lengths are

    1, 4, 7, 14, 19, 32, 39, 58, 67, 92, 103, 134, 147, 184, 199, 242, 259, 308

The last six were found in 2015 with a big optimization solver (integer linear programming). In
2009 David Wilson noticed a pattern and conjectured:

- **odd n:** a(n) = n² − 2n + 4 (so an odd board always misses 2n − 4 squares);
- **even n:** a(n) = n² − n + 2 (an even board misses n − 2 squares).

Nobody had proved either half.

## 4. Why it's hard

To find the best walk by computer you have to rule out every possible walk, and their number
explodes with the board size. That's why the solver stopped at an 18 × 18 board. A formula needs
two different things:

- a **construction**: an actual walk, on every board size, that achieves the length;
- a **proof of impossibility**: an argument that no walk can ever do better.

## 5. What we found

**The odd half of the conjecture is true.** For every odd n, the answer is n² − 2n + 4. So, for
example, **a(19) = 327, a(21) = 403, a(23) = 487**.

**Why nothing can do better (the colouring trick).** Colour each square by whether its column and
row numbers are odd or even. That gives four colours. A horizontal step changes the column's
parity and a vertical step changes the row's, so an alternating walk runs through the four colours
in a fixed cycle. In particular, **every 4 consecutive squares of the walk include a square whose
row and column are both odd**. On an odd board those squares are the scarcest colour. There are
only ((n−1)/2)² of them, and each can be used once, which caps the walk at 4·((n−1)/2)² + 3
squares. That's exactly n² − 2n + 4.

**A walk that achieves it (grow the board, then rotate).** Start with a best walk on the 5 × 5
board. To go from size n to n + 2, add two new columns on the right and two new rows at the
bottom. Extend the walk with a zigzag down the new columns, then a zigzag back along the new rows.
Zigzags are exactly what an alternating walk is good at. This picks up 4n new squares and misses
only 4. Then turn the whole board upside down: the walk now ends in the same kind of spot as
before, ready for the next round. Repeat forever.

**Even boards.** The same trick, starting from the 4 × 4 walk above, gives walks of length
n² − n + 2 on every even board. So the even conjecture can't be *too high*. Whether it can be
beaten is still open.

## 6. How we know it's right

- **A computer check of the construction:** a separate program built the walks for every board
  size up to 301 × 301 and checked every step: alternating, no square repeated, the right length.
- **Agreement with the known terms:** the formulas match all 18 values already in the OEIS,
  including the six found by the optimization solver.
- **A computer-checked proof:** the whole odd case, both the colouring bound and the
  construction, is written in **Lean**. Lean is a proof assistant that accepts a proof only if
  every step is justified down to the axioms of mathematics. It confirms the formula for *every*
  odd n, not just the ones anyone tried.

## 7. Why it matters (honestly)

This is a small puzzle result, not a breakthrough. It still has real value:

- **Half of a 17-year-old conjecture is settled.** Every odd term is now known exactly, forever,
  instead of needing ever-bigger solver runs.
- **A clean example of "impossible" proofs.** The colouring trick is the same style of argument
  as the famous "you can't tile a chessboard with two opposite corners removed using dominoes"
  puzzle. Here it gives an exact answer.
- **The open part is now sharp.** For even boards we know the conjectured value is achievable.
  What's missing is only the impossibility half, a well-defined problem for anyone who likes
  puzzles.
- **It's machine-checked.** Nobody has to trust the argument on faith, including ours.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Self-avoiding walk:** a walk that never visits the same square twice.
- **Colouring argument:** proving something is impossible by colouring squares and counting how
  many of each colour a solution would need.
- **Induction:** proving something for every size by showing how to get from one size to the
  next.
- **Lean:** a proof assistant. It checks every step of a proof down to the axioms of mathematics.
