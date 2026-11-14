# A027851(10) and A029851(10): semigroups of order 10 up to isomorphism, by counting the self-dual ones

**A029851(10) = 613649634393** self-dual (self-converse) semigroups of order 10.

**A027851(10) = 24836001541112970975** semigroups of order 10 up to isomorphism.

Interactive companion page: https://guillochon.github.io/mirror-semigroups/ (every semigroup with up to 5
elements beside its mirror image, and the counting idea on a board you can fill in).

The number of semigroups of order 10 *up to isomorphism or anti-isomorphism* is known:
A001423(10) = 12418001077381302684 (Distler, Jefferson, Kelsey and Kotthoff, CP 2012). The number up to
isomorphism alone was not: it needs the number of semigroups that are isomorphic to their own dual, since

    A027851(n) = 2·A001423(n) − A029851(n).

Both entries stopped at order 9 (Distler and Kelsey 2013). This page computes the self-dual count for
order 10 directly, in about 35 minutes on a desktop PC.

## 1. Burnside over anti-automorphisms

Let S_n × C_2 act on the associative tables on [n], the C_2 by transposing the table. Comparing the orbit
counts of S_n × C_2 and of S_n gives

    A029851(n) = (1/n!) Σ_{π ∈ S_n} A(π),     A(π) = #{associative T on [n] : π(xy) = π(y)π(x) for all x, y}.

A(π) depends only on the cycle type of π. For one representative π of each cycle type (42 for n = 10),

    A(π)/|C(π)| = Σ over the orbits of the centralizer C(π) on those tables of 1/|stabilizer|.

So the whole count is a sum over pairs (semigroup, anti-automorphism) up to isomorphism, each weighted by
1/|its stabilizer|. A semigroup with no anti-automorphism never appears, which is the point: of the
~7·10^14 semigroups of order 10 that are not 3-nilpotent, only ~4.9·10^9 pairs have to be visited.

## 2. The search

`rust/src/main.rs` enumerates the orbits for one π by backtracking over the table:

- **Anti-automorphism.** Filling T[x][y] = v also fills T[πy][πx] = π(v).
- **Associativity.** Each filled cell is propagated through every triple it takes part in
  ((xy)t = x(yt), (tx)y = t(xy), and the two forms where x or y is itself a product).
- **One table per orbit.** The table must be the lexicographically least in its C(π)-orbit (cells read
  diagonal first, then row by row). C(π) is held as an explicit list. Going down the tree, an element whose
  image is already larger is dropped, and a node is abandoned as soon as some image is smaller. At a leaf the
  elements left are exactly the stabilizer.

## 3. Skipping the 3-nilpotent semigroups

Almost all semigroups satisfy S³ = {0}. Their number is a formula (Distler and Mitchell 2012; A383885,
A383886), and so is the number of self-dual ones: 2·A383886(10) − A383885(10) = 608877121317, plus the zero
semigroup. Enumerating 6·10^11 of them would take days, so the search leaves them out:

- If S³ = {0} then x⁴ = 0 for all x. So only diagonals with a unique idempotent z and (x²)² = z for all x
  can contain them. For such a diagonal the search splits in two.
- **z is not a zero.** Then z is the identity of a nontrivial group and the row of z is not constant; this
  is imposed as a constraint.
- **z is a zero.** Then S is nilpotent. The search first chooses the *depth* of every element (the length of
  its longest factorization), requiring some depth ≥ 3, which is exactly S³ ≠ {0}. Depths give strong
  constraints:
  - ab = p ≠ 0 forces depth(p) ≥ depth(a) + depth(b);
  - the depths that occur are exactly 1, …, m;
  - since S^k = S·S^(k−1) = S^(k−1)·S, every element of depth k ≥ 2 must be (generator)·(depth k−1) and
    (depth k−1)·(generator). When no filled cell provides this, the open cells that could are probed.
  The depth vector is part of the object that is made lex-least, and the cells are then read deepest
  products first, so the factorization constraints act early.

Without the last two items the nilpotent branch was 50 times slower.

## 4. Checks

- **Orders 2 to 9, everything.** `selfdual sd n` (no shortcut for S³ = {0}) reproduces A029851(n) for every
  n ≤ 9, and the full rows of A058118 (self-dual semigroups by number of idempotents). At order 9:
  623492664, row 613365656, 8265721, 739317, 410158, 328937, 223226, 113160, 38979, 7510, as in Distler and
  Kelsey's Table 2. The part with S³ = {0} is 606097492 = 2·A383886(9) − A383885(9) + 1.
- **All semigroups.** `selfdual iso n` (same search, no anti-automorphism, group S_n) gives A027851(n) and
  the rows of A058108 for n ≤ 7.
- **The shortcut.** `selfdual sdx n` equals the `sd n` count of semigroups with S³ ≠ {0} for every n ≤ 9
  (17395172 at order 9), by number of idempotents.
- **Commutative semigroups, order 10.** For π = identity the orbits are the commutative semigroups. The run
  found 141655715 with S³ ≠ {0}. Independently, A001426(10) = 3518930337 (Grillet) minus the commutative
  semigroups with S³ = {0}, 3377274622 by a Burnside formula (`check_commutative.py`, which also matches
  orders 8 and 9), is 141655715.
- **Commutative and inverse semigroups through the weights.** Every commutative semigroup and every inverse
  semigroup is self-dual, so the weighted sums Σ count/stab over all 42 classes, restricted to those, must
  be A001426(n) and A001428(n). They are, for n = 7, 8, 9 (17291, 221805, 11545843 and 911, 4637, 26422).
  At order 10 the second run gives 141655715 commutative semigroups with S³ ≠ {0} and 169163 inverse
  semigroups = A001428(10) (Malandro).
- **Two runs.** Order 10 was run twice, with 6 and with 5 threads (different work splits, 35 and 39
  minutes). Every by-idempotents total agrees.
- **Integrality.** For each number of idempotents the sum Σ count/stab is an integer (stabilizers up to
  3628800 occur).
- **What is taken from the literature.** A001423(10) (Distler et al. 2012) and the 3-nilpotent formulas
  (A383885, A383886). A027851(10) is only as good as A001423(10); A029851(10) does not depend on it.

## 5. Results

| Entry | What | New |
|---|---|---|
| A029851 | self-converse semigroups | a(10) = 613649634393 |
| A027851 | semigroups up to isomorphism | a(10) = 24836001541112970975 |
| A058118 | self-converse semigroups of order n with k idempotents | row 10: 611179527739, 2419403101, 36082031, 5455752, 3605531, 2735993, 1712924, 806327, 258856, 46139 |
| A058119 | self-converse, 1 idempotent | a(8)–a(10) = 1937266, 613365656, 611179527739 |
| A058120 | self-converse, 2 idempotents | a(8)–a(10) = 123988, 8265721, 2419403101 |
| A058121 | self-converse, 3 idempotents | a(8)–a(10) = 49969, 739317, 36082031 |
| A058122 | self-converse bands | a(8)–a(10) = 1376, 7510, 46139 |
| A058109 | semigroups with 1 idempotent, up to isomorphism | a(8)–a(10) = 3667785000, 105952488687468, 24835238539006265743 |
| A058110 | 2 idempotents | a(8)–a(10) = 9207430, 25412267163, 762557534616451 |
| A058111 | 3 idempotents | a(8)–a(10) = 1150085, 136799017, 439138761619 |
| A058112 | bands up to isomorphism | a(10) = 14020041 |

The order 8 and 9 values of the column entries are already in the triangles A058118 and A058108 (and are
reproduced here); the order 10 values are new. The "up to isomorphism" values are
2·(up to equivalence) − (self-dual), with A002786(10), A002787(10), A005591(10) and A002788(10).

Of the 613649634393 self-dual semigroups of order 10, 608877121318 have S³ = {0} and 4772513075 do not.
The search visited 4883233097 (semigroup, anti-automorphism) classes in 1.15·10^10 nodes; the four
largest cycle types (the identity and one, two or three transpositions) hold 99.9% of them.

## 6. Reproduce

    cd rust && cargo build --release
    target/release/selfdual sd 9 6 > out9.txt        # everything, order 9: 3 min
    target/release/selfdual sdx 9 6 > outx9.txt      # without S^3 = 0: 18 s
    target/release/selfdual sdx 10 6 > outx_n10.txt 2> logx_n10.txt
    python combine.py out9.txt
    python derive.py outx_n10.txt
    python check_commutative.py logx_n10.txt

Output lines are `class <i> <cycle type> k=<idempotents> nil3= com= inv= stab=<s> count=<c>`; a count of
semigroups is Σ count/stab.

## Files

- `rust/src/main.rs`: the search (`sd`: everything; `sdx`: without S³ = {0}; `iso`: all semigroups up to
  isomorphism, a check for n ≤ 7).
- `combine.py`: sums an output file. `derive.py`: the OEIS values. `check_commutative.py`: Burnside count
  of commutative semigroups with S³ = 0, and the identity-class check.
- `out_n*.txt`, `outx_n*.txt`, `log*.txt`: outputs and logs of the runs.
