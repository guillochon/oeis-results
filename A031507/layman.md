# A031507 and A031508 in plain English

*What these sequences count, why anyone cares, and what we found.*

---

## 1. The idea

Take the equation **y² = x³ + k** for a whole number k. Its graph is a smooth curve, a so-called *elliptic curve*.
The question that has fascinated mathematicians since antiquity is: which points on it have **rational coordinates**
(fractions)?

Here is a surprise. You can often build new rational points from old ones: draw the line through two of them, and
it meets the curve in a third rational point. All the rational points of the curve can be generated this way from a
finite starting set. The **rank** of the curve is the number of truly independent starting points you need, not
counting a few "torsion" points of finite order.

- For y² = x³ + 1 there are only finitely many rational points (rank 0).
- For y² = x³ + 2 the point (−1, 1) generates infinitely many (rank 1).
- For y² = x³ + 15 you need two independent points, (1, 4) and (109, 1138) (rank 2).

**A031507(n)** is the smallest k > 0 for which y² = x³ + k has rank exactly n. **A031508(n)** is the same for
y² = x³ − k. **A373795** takes the smaller of the two.

## 2. Where this shows up

- **Mordell's equation.** y² = x³ + k is one of the oldest Diophantine equations. Its rational points were studied by
  Fermat, and Mordell proved the "finitely generated" theorem for all elliptic curves in 1922.
- **The rank problem.** How large can the rank of an elliptic curve be? Nobody knows. Mathematicians hunt for record
  ranks, and these sequences ask the cleanest version: how soon does each rank first appear?
- **A famous open problem.** The Birch and Swinnerton-Dyer conjecture, one of the million-dollar Millennium
  Prize Problems, is about exactly this rank. Many rank computations quietly assume it (or the Riemann hypothesis for
  number fields, "GRH"). Ours does not.

## 3. What was known

The OEIS listed ranks 0 to 6:

- A031507: 1, 2, 15, 113, 2089, 66265, 1358556
- A031508: 1, 2, 11, 174, 2351, 28279, 975379

The rank-6 values were proved in 2001 by Tom Womack, who ran the program mwrank on every curve up to about
1.4 million. For rank 7 he found the curves **y² = x³ + 47550317** and **y² = x³ − 56877643**. But he could not
prove that no smaller k works. He estimated that checking the roughly 100 million smaller curves would take
**"one year running on 500 1GHz computers"**. So a(7) stayed a conjecture for 24 years.

## 4. Why it's hard

Computing the rank of one curve is itself hard. The standard method ("2-descent") needs detailed information about
a number field attached to the curve (its *class group*). Computers can only get that information quickly by
assuming GRH; proving it without GRH is slow. Doing this 100 million times is out of reach.

## 5. What we found

**A031507(7) = 47550317 and A031508(7) = 56877643**, so also A373795(7) = 47550317. Both are **proved with no
unproved conjectures**.

The trick is to **rule curves out cheaply, in stages**, keeping the expensive method for a handful of stubborn cases:

1. **A cheap ceiling on the rank.** These curves have a special symmetry (a "3-isogeny") that gives a quick upper
   bound on the rank. The bound depends on the class groups of two *quadratic* fields, Q(√k) and Q(√−3k), which are
   much simpler than the cubic field the standard method needs. Better still, we only need the "3-part" of those
   class groups. By a classical theorem of Hasse, that can be read off from **how many cubic number fields have a
   given discriminant**. So we counted all cubic fields up to discriminant 683 million, and that gave the needed
   information for every k at once, without assuming anything. With that table, a compiled program checked all
   114 million curves in **2.5 seconds**. All but about 72,000 were ruled out.
2. **A balance law.** A formula of Cassels (1965) relates the two halves of the bound exactly. That rules out about
   40,000 more.
3. **Explicit counterexamples.** For each remaining curve we looked for a concrete number that the bound "counts"
   but that fails a simple test modulo some prime. Each one found lowers the ceiling by 2. That leaves 264 curves.
4. **The expensive method, made rigorous.** For those 264 curves we ran the standard 2-descent, but first *proved* the
   class-group information it relies on, using PARI's `bnfcertify`. This took about 4 CPU-hours; one curve alone
   took an hour. Every one of them has rank at most 5.

Finally, the two rank-7 curves really have rank 7. We have 7 explicit independent points on each, and the cheap
ceiling from step 1 is exactly 7 for both.

The whole computation took about 5 CPU-hours on a home PC, against Womack's estimate of 500 machine-years for his
approach.

## 6. How we know it's right

- **It reproduces the known answer.** The same pipeline, aimed at rank 6, re-proves Womack's a(6) for both
  sequences in a few seconds, by a completely different method.
- **Every ceiling was tested against reality.** On thousands of small curves whose rank is known, no ceiling ever came
  out below the true rank. As a control, we deliberately introduced a sign error; it was caught thousands of times.
- **The "counterexample" test was tested too.** It never rejected a number coming from an actual rational point
  (427 checked).
- **The class-group table was cross-checked.** It agrees with PARI's standard (GRH-based) class group computation for
  all 1.2 million discriminants up to 2 million. Every one of its 120 million entries also has the exact form
  Hasse's theorem demands.
- **One small caveat.** Like almost all computer algebra, the factoring uses a primality test (BPSW) that is proven
  for numbers below 2⁶⁴ and has no known failure beyond.

## 7. Why it matters (honestly)

This is a solid, citable result. It is not a breakthrough in the theory of ranks.

- **It closes a 24-year-old gap** in two well-known OEIS entries by settling Womack's conjectured values.
- **It is unconditional.** Many results in this area assume BSD or GRH. For example, the analogous minimality
  results for the related curves x³ + y³ = k are known only conditionally. This one assumes neither.
- **It is cheap and reusable.** The method (a symmetry-based ceiling, a cubic-field table, and expensive rigorous
  work only at the end) is far cheaper than running the full method on every curve. It should extend to rank 8, whose
  candidates are about 30 times larger.
- **Honest limits.** The rank-7 curves themselves were already found by Womack. Our contribution is the proof that
  nothing smaller works.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number sequences.
- **Elliptic curve:** a curve like y² = x³ + k; its rational points form a group.
- **Rank:** the number of independent rational points needed to generate all of them (apart from finitely many).
- **Class group:** a finite group measuring how badly unique factorization fails in a number field.
- **GRH / BSD:** two famous unproved conjectures (the generalized Riemann hypothesis and the Birch and Swinnerton-Dyer
  conjecture), often assumed in computations like this one.
