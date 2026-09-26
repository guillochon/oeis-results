# A236553 in plain English

*What this sequence counts, what we found, and how we know it's right.*

---

## 1. The question

**Split quaternions** are numbers of the form a + b·i + c·j + d·k with i² = j² = +1. They're a
cousin of Hamilton's quaternions, which describe 3D rotations, and are really 2×2 matrices in
disguise. Do the arithmetic **modulo n**, like clock arithmetic with an n-hour clock, so that a,
b, c, d are remainders after dividing by n.

An **involution** is an X with **X² = 1**: something that undoes itself when done twice, like a
mirror reflection. **A236553(n)** counts the involutions when arithmetic is done modulo n.
Equivalently, it counts the solutions of

    a² − b² + c² + d² = 1,   2ab = 0,   2ac = 0,   2ad = 0      (modulo n)

The OEIS listed the first 49 values (1, 8, 14, 64, 32, 112, 58, 288, …), found by a computer
search, but no formula.

## 2. The first simplification: build from prime powers

Clock arithmetic has a famous property, the **Chinese remainder theorem**. Arithmetic modulo 12
is "the same as" doing arithmetic modulo 4 and modulo 3 side by side, because 4 and 3 share no
factor. So a solution modulo 12 is the same thing as a solution modulo 4 paired with a solution
modulo 3, and

    A236553(12) = A236553(4) × A236553(3) = 64 × 14 = 896.

Sequences with this property are called *multiplicative*, and the OEIS already tagged this one
`mult`. It means we only need a formula for **prime powers**: 3, 9, 27, …, 5, 25, …, and the
powers of 2.

## 3. What we found

> For an odd prime p: **A236553(pᵏ) = p^(2k−1) · (p + 1) + 2.**
>
> For powers of 2: **8, 64, and then 2^(2k+2) + 32** from 8 on (this is sequence A236554).

For example, A236553(5) = 5 × 6 + 2 = 32 and A236553(9) = 27 × 4 + 2 = 110, both matching the
OEIS. Combined with multiplicativity, this gives **every** term. The ten-thousandth term takes a
fraction of a second instead of a huge search.

## 4. Why, roughly

For an odd prime power, split the solutions by the first coordinate a.

- **a is "invertible"** (not a multiple of p): the other equations force b = c = d = 0 and
  a = ±1. That's the **"+2"**.
- **a is a nonzero multiple of p:** impossible. Reducing everything modulo p gives 0 = 1.
- **a = 0:** we have to count solutions of c² + d² − b² = 1.
  - **Modulo p itself** there are exactly **p² + p**, found by a change of variables that turns
    the equation into c² + u·w = 1, which is easy to count.
  - **Going up from pʲ to pʲ⁺¹**, each solution has p³ possible "children", and exactly p² of
    them are solutions. (This is a counted version of *Hensel's lemma*, a classical tool for
    lifting solutions to higher powers.) So the count gets multiplied by p² at each step, giving
    p^(2k−1)(p + 1).

## 5. How we know it's right

1. **Brute force:** a direct count from the definition, for every n up to 256, agrees with the
   formula. So do all 49 values in the OEIS.
2. **A written proof**, in `README.md`.
3. **A computer-checked proof in Lean.** Lean is a proof assistant, a programming language where
   every logical step has to be justified down to the axioms of mathematics. It verifies the
   formula for **all** n, not just the ones we tested. It rests only on the three standard axioms
   that all of Lean's mathematics library uses.

A sibling sequence, **A227867**, asks the same question for ordinary quaternions (i² = j² = −1).
Surprisingly, it has exactly the same values at odd n. The two differ only at powers of 2, and
there's a nice reason why (see its own explanation).
