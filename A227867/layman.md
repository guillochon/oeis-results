# A227867 in plain English

*What this sequence counts, what we found, and how we know it's right.*

---

## 1. The question

**Quaternions** were invented by Hamilton in 1843. They are numbers a + b·i + c·j + d·k with
i² = j² = k² = −1, and they are how computers describe 3D rotations. A **Lipschitz quaternion** is
one whose coefficients a, b, c, d are whole numbers.

Now do the arithmetic **modulo n**, like clock arithmetic with an n-hour clock, and ask: how many
quaternions X satisfy **X² = 1**? These are the "involutions", things that undo themselves when
done twice, like a mirror reflection. That count is **A227867(n)**. Written out, it counts the
solutions of

    a² − b² − c² − d² = 1,   2ab = 0,   2ac = 0,   2ad = 0      (modulo n)

The OEIS had 52 values from computer searches (1, 8, 14, 32, 32, 112, 58, 32, …), but no formula.

## 2. What we found

The sequence is **multiplicative**: the value at 12 is the value at 4 times the value at 3, and
so on (Chinese remainder theorem). So everything comes down to prime powers:

> For an odd prime p: **A227867(pᵏ) = p^(2k−1) · (p + 1) + 2.**
>
> For powers of 2: **A227867(2) = 8, and A227867(2ᵏ) = 32 for every k ≥ 2.**

Two things stand out:

- **At odd numbers it matches A236553 exactly.** A236553 is the same question for the "split"
  quaternions, where i² = j² = +1. At odd moduli, the sign of i² makes no difference to the
  count.
- **At powers of 2 it freezes at 32.** Modulo 4, 8, 16, 1024, or 2¹⁰⁰⁰, there are always exactly
  32 solutions, while the split version keeps growing (288, 1056, 4128, …).

## 3. Why

**Odd primes.** The count comes down to the number of points on a "sphere", b² + c² + d² = −1,
modulo p. The trick is that −1 can always be written as a sum of two squares modulo p, say
−1 = s² + t². (For p = 7: −1 ≡ 6, and 2² + 3² = 13 ≡ 6 modulo 7.) Using s and t you can rotate the coordinates so that the sphere turns into
exactly the surface counted for A236553. So both have p² + p points modulo p, and from there the
counts grow in the same way.

**Powers of 2: the "7 mod 8" rule.** Squares modulo 8 can only be 0, 1 or 4. Add up three of
them and you can get 0 through 6, but **never 7**. This fact goes back to Legendre and Gauss:
numbers of the form 8m + 7 are never a sum of three squares. For A227867, half of the potential
solutions would need b² + c² + d² ≡ 7 (mod 8), so they simply don't exist. What's left is
always 4 × 8 = 32 solutions, whatever the power of 2.

## 4. How we know it's right

1. **Brute force:** a direct count from the definition agrees with the formula for every n up
   to 256, and with all 52 OEIS values.
2. **A written proof**, in `README.md`.
3. **A computer-checked proof in Lean**, a proof assistant that accepts an argument only if every
   step is justified down to the axioms. It confirms the formula for **all** n. The "never 7
   mod 8" fact is checked by having Lean try all 512 combinations modulo 8.

## 5. Why it matters (honestly)

This is a modest result, not a breakthrough. It still has real value:

- **A complete answer.** The OEIS had 52 values from computer searches and no formula. Now every
  term follows from a short product over the prime factors of n.
- **It explains why two sequences agree.** At odd n, A227867 equals A236553, the same count for
  the "split" quaternions. Both equal the number of 2×2 matrices X with X² = I modulo n. The
  reason is that, modulo an odd number, ordinary quaternions and split quaternions are both
  secretly the same ring of 2×2 matrices. That's a real structural fact, now visible in the
  numbers.
- **An old theorem shows up.** The freeze at 32 for powers of 2 is the three-square theorem of
  Legendre and Gauss ("8m + 7 is never a sum of three squares") appearing in a new place.
- **It's machine-checked.** The proof is verified in Lean, down to the axioms, including the
  "never 7 mod 8" fact.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Quaternion:** a number a + b·i + c·j + d·k with i² = j² = k² = −1. Lipschitz quaternions have
  whole-number coefficients.
- **Modulo n:** clock arithmetic, keeping only the remainder after dividing by n.
- **Multiplicative:** the value at a product of coprime numbers is the product of the values.
- **Lean:** a proof assistant. It checks every step of a proof down to the axioms of mathematics.
