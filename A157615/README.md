# A157615 — Longest alternating (H/V) self-avoiding path on an n × n board

**Status (2026-09-26):**
- **Odd n: solved, and formally verified in Lean 4 + Mathlib** (`lean/OeisLean/A157615.lean`, theorem `A157615.a_odd`; axioms: propext, Classical.choice, Quot.sound only). a(n) = n² − 2n + 4 for every odd n ≥ 3. The four-colouring argument gives the upper bound (below), and an explicit frame-induction construction (`construct.py`) gives a path of that length for every odd n.
- **Even n:** the construction gives a path of length n² − n + 2 for every even n ≥ 4, so a(n) ≥ n² − n + 2, the conjectured value. **The matching upper bound is still open.** It is known only by ILP for n ≤ 18.
**Known terms:** a(1)..a(18) = 1, 4, 7, 14, 19, 32, 39, 58, 67, 92, 103, 134, 147, 184, 199, 242, 259, 308. a(13)..a(18) come from integer linear programming (Rob Pratt, 2015).
**Conjecture (David Wilson, in the entry):** a(n) = n² − n + 2 for even n; a(n) = n² − 2n + 4 for odd n > 1.
**Targets:**
1. Prove the conjecture. That settles every term and lets `hard`/`more` come off.
2. As a near-term contribution: a(19) = 327 and a(21) = 403 need only a construction, since the odd upper bound is already proved. a(20) = 382 also needs the even bound or a solver proof.

**Bonus siblings (also `hard,more`), with the same moves and the same colouring argument:**
- **A157616**: the path must start at a corner. Its odd-n terms are exactly one less than A157615's.
- **A157617**: closed cycles instead of paths.

## What is being counted
You walk on the squares of an n × n board. Each step goes to a horizontally or vertically adjacent square, and **the steps must alternate**: horizontal, vertical, horizontal, … You may never revisit a square. a(n) is the largest number of squares such a path can cover.

Because the steps alternate, **the path turns 90° at every square**. It is a zigzag or staircase that never goes straight. `brute.py` searches exhaustively and reproduces a(1..6). It prints optimal paths, for example n = 4 (14 of 16 squares, missing two corners):

    S # # E
    # # # #
    # # # #
    . # # .

## Proved: the odd-n upper bound (four-colouring)
Colour each square by (x mod 2, y mod 2), which gives 4 colours. A horizontal step flips the x-parity and a vertical step flips the y-parity. So an alternating path cycles through the colours in a fixed order, period 4:

    (0,0) → (1,0) → (1,1) → (0,1) → (0,0) → …

Every colour therefore appears ⌊L/4⌋ or ⌈L/4⌉ times in a path of length L, so L ≤ 4·m + 3, where m is the size of the smallest colour class.

For **odd n**, the smallest class is the squares with both coordinates odd: m = ((n−1)/2)². That gives

    L ≤ (n−1)² + 3 = n² − 2n + 4,

**exactly the conjectured value.** So for odd n only a construction is missing.

For **even n**, all four classes have n²/4 squares, so the colouring gives only the useless bound L ≤ n² + 3. The conjectured deficit of n − 2 squares must come from somewhere else. That is the heart of the problem.

## The construction (proves all odd terms)
`construct.py` builds the paths explicitly and validates them with `path.py` for every n ≤ 301.

**Odd n.** Invariant: P(n) has length n² − 2n + 4 and ends (n−1, 0) → (n−1, 1), on the right edge, arriving downward. To go from n to n + 2, add a width-2 L-shaped frame (columns n, n+1 and rows n, n+1) and append the route:
1. (n, 1);
2. zigzag down the right strip: in rows r = 2..n, visit (n, r), (n+1, r) if r is even, and (n+1, r), (n, r) if r is odd;
3. (n, n+1);
4. zigzag left along the bottom strip: in columns c = n−1..0, visit (c, n+1), (c, n) if n−1−c is even, and (c, n), (c, n+1) otherwise.

The route adds exactly 4n squares, missing only (n, 0), (n+1, 0), (n+1, 1) and (n+1, n+1). It ends (0, n+1) → (0, n), and rotating the new board 180° restores the invariant.
- Base case P(5): 19 squares, found by exhaustive search.
- Each join alternates correctly: the vertical arrival at (n−1, 1), the turn from the right strip into the bottom strip, and the start of each zigzag.
- The route stays inside the new frame, so it is self-avoiding.

Together with the colouring bound, this gives **a(n) = n² − 2n + 4 for all odd n ≥ 3**. The first new terms are **a(19) = 327, a(21) = 403, a(23) = 487**.

A consequence of the colouring argument: an optimal odd path visits every (odd, odd) square, and starts and ends on "mixed" squares (one coordinate odd, the other even). This is why the invariant endpoint (n−1, 1) is a mixed square.

**Even n.** Invariant: E(n) has length n² − n + 2 and ends (n−1, 1) → (n−1, 0), at the top-right corner, arriving upward. The frame route has the same shape: (n, 0); then (n, r), (n+1, r) for odd r and (n+1, r), (n, r) for even r, over r = 1..n; then (n, n+1); then the same bottom zigzag. It adds 4n + 2 squares, missing only (n+1, 0) and (n+1, n+1). The base case is E(4), 14 squares. So **a(n) ≥ n² − n + 2 for all even n ≥ 4.**

## Still open: the even upper bound
The colouring gives nothing for even n, because all four classes have size n²/4. A reformulation that might help: in any valid path, each interior square has exactly one horizontal and one vertical edge. So the visited set splits into horizontal dominoes within rows, and also into vertical dominoes within columns (up to the two endpoints). The union of the two matchings must be a single path, with no alternating cycles. For even n the whole board admits both matchings, but their union is then a set of 4-cycles (the 2×2 blocks). The deficit of n − 2 squares must be the cost of breaking every alternating cycle.

## Clever ideas to pursue
1. **Done:** odd constructions by frame induction (see above). *Original idea:* Going from n to n + 2 must add exactly 4n squares, one frame of width 2 minus 4. Get optimal paths for n = 5, 7, 9, 11 from a solver, look for a rule that splices a width-2 frame into an existing path, and prove the rule preserves alternation and self-avoidance. This gives every odd term. **a(19) = 327 is then immediately provable**: exhibit a length-327 path and cite the colouring bound.
2. **Done:** even constructions by the same method (n → n + 2 adds 4n + 2).
3. **Boundary rigidity lemma (a lead for the even upper bound).** Look at the top row. Each square there that is interior to the path uses its only vertical neighbour (below) and one horizontal neighbour. So the visited top-row squares pair up into horizontal dominoes, each with both squares stepping down. Below a domino {x, x+1}, the square (x,1) cannot connect to (x+1,1), because that closes a 4-cycle. So (x,1) must go left and (x+1,1) must go right. A domino starting at a corner therefore forces a path endpoint, and the path has only two. This rigidity spreads diagonally inward, which fits a linear (n − 2) deficit. Next step: turn it into a counting argument, for example with a discharging scheme that charges each missed square to a boundary obstruction.
4. **Let the computer suggest the proof (LP duality).** Write the problem as an integer program (as Pratt did). Add valid local constraints: alternation, degree ≤ 2, the colour-class balance, the domino rule from idea 3. Solve the **LP relaxation** for even n = 4…14 and inspect the **dual solution**. If a structured family of dual weights certifies n² − n + 2, those weights are a proof skeleton: "every square gets weight w, every valid path has total weight ≤ bound". This is a strong way for AI and computation to find a human-readable proof.
5. **A finer colouring for even n.** Look for a colouring or weighting mod 4 (or along diagonals x ± y) in which the path's cyclic structure forces an imbalance of about n/2 per pair of sides. Search small weightings by computer against the exact values for n = 4, 6, 8.
6. **Generalize to m × n rectangles.** For m odd and n even the colouring bound is L ≤ (m−1)·n + 3; for both odd it is (m−1)(n−1) + 3. Exact values for small rectangles (from `brute.py`, adapted) could reveal the true formula, and a rectangle statement is often easier to prove by induction than the square one.

## Files
| File | What it is |
|---|---|
| `README.md` | This write-up: the proofs, the open even case, and ideas. |
| `layman.md` | Plain-English explanation. |
| `path.py` | Independent validator for alternating self-avoiding paths. |
| `construct.py` | The frame-induction constructions for odd and even n, validated for all n ≤ 301. |
| `brute.py` | Exhaustive search for n ≤ 6 (reproduces a(1..6)). |
| `search.py` | Randomized Warnsdorff search; reaches the conjectured length for n ≤ 12. |
| `lean/` | Lean 4 + Mathlib formal proof of the odd case (`OeisLean/A157615.lean`). Build with `cd lean && lake exe cache get && lake build`. |

## Tools needed
- A SAT/CP solver: `pip install ortools` (CP-SAT) or `python-sat`. Neither is installed yet.
- Model: Boolean x[square] for "on the path", Boolean edge variables e[h/v edge], degree ≤ 2 at every square, the alternation constraint (at a degree-2 square, exactly one horizontal and one vertical edge), and exactly two degree-1 squares. Connectivity or no-subcycles is the tricky part: add lazy cycle cuts, or use CP-SAT's `AddCircuit` with a dummy node.

## Verification plan
- Any constructed path is checked by a short independent validator: alternation, self-avoidance, length.
- The even upper bound, once found, is checked against solver optima for n ≤ 18 (all known).
- Extend `brute.py` to n = 7 with better pruning (the colour-count bound) to cross-check the solver model.

## Log
- 2026-09-26: `brute.py` confirms a(1..6). Found that the four-colouring argument proves the odd-n upper bound n² − 2n + 4 exactly, so the odd case reduces to constructions.
- 2026-09-26: Odd case formalized in Lean (`a_odd`: the colouring upper bound via disjoint 4-windows, the frame route as an explicit piecewise function, the induction via extend-then-rotate, and explicit 3×3 and 5×5 base cases).
- 2026-09-26: `search.py` (randomized Warnsdorff DFS) reaches the conjectured length for every n ≤ 12 within seconds. A search for a self-similar frame step found one from 5×5 (odd) and from 4×4 (even). Generalized in `construct.py`, validated for all n ≤ 301. **The odd case is solved**, and the even lower bound is proved.
