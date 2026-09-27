# A109457, A109458, A109459 in plain English

*What these sequences count, why anyone cares, and what we found.*

---

## 1. The idea

A logic puzzle often comes as a list of simple rules such as "if the light is on then the door is
open" or "at least one of Alice and Bob must attend". Each rule mentions at most two things. A
collection of such two-item rules is called a **2-SAT instance** (or a **Krom formula**, after
Melven Krom, who studied them in 1967). 2-SAT is famous because it is *easy*: a computer can
decide in linear time whether all the rules can be satisfied at once, while the same question
with three-item rules (3-SAT) is the textbook hard problem.

## 2. Functions, not formulas

Many different rule lists say the same thing: "if A then B" together with "if B then A" is the
same as "A and B are equal". What matters is the set of situations the rules allow. That set is
a **Boolean function** of the n variables, and a function that can be described by two-item
rules is a **Krom function**. There is a neat way to recognise one without any formula: take any
three allowed situations and, for each variable, let the majority of the three decide. A function
is Krom exactly when this "majority vote" of allowed situations is always allowed again.

## 3. What the three sequences count

- **A109457(n)** counts the Krom functions of n variables. For n = 2 there are 16 (every function
  of two variables is Krom).
- **A109458(n)** counts them up to *renaming* the variables.
- **A109459(n)** counts them up to renaming *and* flipping variables (swapping "on" and "off").

These are the numbers that tell you how many genuinely different 2-SAT problems there are on n
variables. Donald Knuth computed them in 2005 for his book *The Art of Computer Programming*
(Volume 4A, and Problem 191 of the *Satisfiability* fascicle): A109457 and A109458 to n = 8,
A109459 to n = 7. They were left marked "hard, more" in the OEIS.

## 4. Why the direct approach stalls

The labelled count grows fast: about 2·10^12 functions for n = 8. Listing them one by one, even
at a billion per second, was already out of reach for Knuth at n = 9 and is hopeless at n = 10.

## 5. The trick: describe a Krom function by its skeleton

A satisfiable 2-SAT instance, once you add every rule it implies, has a rigid shape:

1. some variables are simply **forced** to a fixed value;
2. some groups of variables are forced to move **together** (equal or opposite);
3. what remains is a **web of one-way implications** among the groups, with a built-in mirror
   symmetry: "if A then B" always comes with "if not-B then not-A".

We call the third layer a *skew poset*: a partial order on the 2r "signed" groups that is turned
upside down by negation, and in which no group implies its own opposite. Once you see this, the
counting problem changes character. Instead of listing functions you list skeletons, and there
are very few of them: 8.8 million for r = 9, versus 1.5·10^15 labelled functions.

For each skeleton we ask a computer program to find all its symmetries (using Brendan McKay's
`nauty` package, the standard tool for such questions). Classic counting formulas of Burnside,
Pólya and Euler then turn the list of skeletons and their symmetries into all three sequences at
once: the labelled count, the count up to renaming, and the count up to renaming and flipping.

## 6. What we found

| n | A109457 (labelled) | A109458 (renaming) | A109459 (renaming + flipping) |
|---|---|---|---|
| 8 | 2061662323954 (Knuth) | 62154403 (Knuth) | **315080 (new)** |
| 9 | **1517417860795700 (new)** | **4787612373 (new)** | **11073674 (new)** |

The whole computation up to n = 9 takes 20 seconds on a desktop PC.

## 7. How we know it is right

- The same program reproduces all 26 values Knuth published, including the 13-digit labelled
  count for n = 8, which depends on every one of the 227348 skeletons having exactly the right
  symmetry group.
- The two "up to symmetry" sequences were recomputed by a second, deliberately different
  program that lists the coloured skeletons directly and never touches the Pólya machinery. The
  two programs agree on every number they both produce.
- Several internal consistency checks (group orders, divisibility, an exponential-formula
  identity between connected and all skeletons) are asserted at run time.

## 8. Where to look

- The OEIS entries: [A109457](https://oeis.org/A109457), [A109458](https://oeis.org/A109458),
  [A109459](https://oeis.org/A109459).
- M. R. Krom, *The decision problem for a class of first-order formulas in which all
  disjunctions are binary*, 1967.
- D. E. Knuth, *The Art of Computer Programming*, Vol. 4A, Section 7.1.1, and Vol. 4 Fascicle 6
  (*Satisfiability*), Problem 191.
- B. D. McKay and A. Piperno, *Practical graph isomorphism, II* (the nauty package).
