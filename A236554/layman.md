# A236554 in plain English

*What this sequence counts, what we found, and how we know it's right.*

---

## 1. Clock arithmetic

On a clock, 9 o'clock plus 5 hours is 2 o'clock, not 14. Mathematicians call this arithmetic
**modulo 12**: every number is replaced by its remainder after dividing by 12. You can add,
subtract and multiply this way, and the usual rules still work. You can do the same with any
modulus. This sequence uses the powers of two: 2, 4, 8, 16, 32, …

## 2. Quaternions, and a "split" cousin

Ordinary numbers live on a line. Complex numbers add one extra direction, *i*, with i² = −1.
**Quaternions**, invented by Hamilton in 1843, add three extra directions: *i*, *j* and *k*. A
quaternion looks like

    a + b·i + c·j + d·k

and multiplication follows a small table of rules (i·j = k, j·i = −k, and so on). Quaternions are
used every day to describe **rotations in 3D**, in video games, spacecraft attitude control and
phone orientation sensors.

This sequence uses a close cousin, the **split quaternions**, where i² = +1 and j² = +1 instead of
−1. They show up in the geometry of spacetime (special relativity) and in studying 2×2 matrices.
In fact, split quaternions are really just 2×2 matrices in disguise.

We do split-quaternion arithmetic **modulo 2ⁿ**: every coefficient a, b, c, d is a remainder
after dividing by 2ⁿ.

## 3. What A236554 counts

An **involution** is something that undoes itself when done twice. Flipping a light switch, or
reflecting in a mirror, are examples. In algebra it means an element X with **X² = 1**. Among
ordinary numbers only +1 and −1 qualify. In a quaternion ring there can be many more.

**A236554(n)** is the number of split quaternions X, with coefficients modulo 2ⁿ, that satisfy
X² = 1. Writing out X² = 1 gives four equations in a, b, c, d:

    a² − b² + c² + d² = 1,   2ab = 0,   2ac = 0,   2ad = 0      (all modulo 2ⁿ)

## 4. What was known, and what we found

The OEIS listed 8 terms, computed in 2014 by a program that tried every possible (a, b, c, d):

    8, 64, 288, 1056, 4128, 16416, 65568, 262176

The entry was tagged `hard` ("terms are hard to compute"), because each new term means checking
16 times as many candidates as the one before.

It turns out the answer has a simple formula. From n = 3 on:

> **A236554(n) = 2^(2n+2) + 32**, that is, 4^(n+1) + 32.

For example, 4⁴ + 32 = 288 and 4⁹ + 32 = 262176. So the next terms are 1048608, 4194336,
16777248, and so on for every n. The sequence isn't hard after all. The difficulty was an artifact
of the brute-force program.

## 5. Why the formula is true, roughly

Split the solutions by whether **a is odd or even**.

- **a odd:** the equations 2ab = 0 and so on pin b, c and d down to just 2 choices each, and a
  itself turns out to have exactly 4 choices. That gives 4 × 2 × 2 × 2 = **32** solutions, the
  "+32", no matter how big n is.
- **a even:** almost every even a leads to a contradiction. Only a = 0 and one other value
  survive, and then you have to count the solutions of c² + d² − b² = 1 modulo 2ⁿ. The key step
  shows that **going from 2ⁿ to 2ⁿ⁺¹ multiplies that count by exactly 4**. Every solution has
  8 "children" one level up, and a clever pairing shows that exactly half of them work. Starting
  from the count at n = 3, this gives the 4^(n+1) part.

## 6. How we know it's right

We checked three independent ways:

1. **Brute force.** A direct count from the definition reproduces all 8 known terms and five new
   ones. That count uses none of the proof's ideas.
2. **A written proof**, in `README.md`.
3. **A computer-checked proof.** The whole argument is written in **Lean**, a "proof assistant":
   a programming language where every logical step must be justified, down to the axioms of
   mathematics. Lean refuses to accept the proof unless every step is correct. It confirms the
   formula for *every* n ≥ 3, not just the ones we tried. So the result doesn't depend on anyone
   (human or AI) having checked the argument carefully.

This sequence also opened the door to two relatives, A236553 and A227867, the same question with
any modulus. Both now have proved formulas too.
