# A236553 — Involutions in the quaternion ring over Z/nZ with i^2 = j^2 = 1: a closed form

**Status:** **proved and formally verified in Lean 4 + Mathlib** (`lean/OeisLean/A236553.lean`, theorem `A236553.A_closed_form`; axioms: propext, Classical.choice, Quot.sound only). Checked independently by direct count for every n ≤ 256 and against all OEIS data.

**Sequence:** A236553(n) is the number of solutions mod n of

    a² − b² + c² + d² = 1,   2ab = 2ac = 2ad = 0.

These are the involutions X² = 1 in the quaternion ring over Z/n with i² = j² = +1. Keyword `mult`; the entry has no formula.

**Related:** A236554(n) = A236553(2^n) (see `../A236554/`). A227867, the Lipschitz quaternions i² = j² = −1, has the same odd part and a different 2-part (see `../A227867/`).

## Result
A236553 is multiplicative, with

| | odd p, k ≥ 1 | 2 | 4 | 2^k, k ≥ 3 |
|---|---|---|---|---|
| **A236553(p^k)** | p^(2k−1)·(p+1) + 2 | 8 | 64 | 2^(2k+2) + 32 |

For odd n this reads

    a(n) = Π_{p^k || n} ( p^(2k−1)·(p+1) + 2 ).

Here p^(2k−1)(p+1) = p^(2k)(1 + 1/p), which is the p-part of n·ψ(n) (ψ is the Dedekind psi function, A001615). So a(n) is "n·ψ(n) with +2 at each prime". For squarefree odd n, a(n) = Π_{p|n} (p² + p + 2).

Spot checks: a(3) = 3·4 + 2 = 14, a(9) = 27·4 + 2 = 110, a(12) = 64·14 = 896. All match the OEIS data.

## Proof
Write N = p^k and Q(b, c, d) = c² + d² − b². The system is a² + Q(b, c, d) ≡ 1 with 2ab ≡ 2ac ≡ 2ad ≡ 0.

**Multiplicativity.** By the Chinese remainder theorem, Z/mn ≅ Z/m × Z/n for coprime m, n. A tuple satisfies the equations mod mn exactly when it satisfies them mod m and mod n.

### Odd p
Since 2 is invertible, the side conditions are ab ≡ ac ≡ ad ≡ 0.
- **a a unit:** b = c = d = 0 and a² ≡ 1. That gives **2** solutions (±1).
- **a a non-unit, a ≠ 0:** if b were a unit, ab ≡ 0 would force a = 0. So b, and likewise c and d, are divisible by p. Then all of a, b, c, d are ≡ 0 (mod p), and the equation says 0 ≡ 1 (mod p). **No** solutions.
- **a = 0:** b, c, d are free, and we need M = #{Q(b, c, d) ≡ 1 mod p^k}.

**Counting mod p.** Set u = d − b and w = d + b, which is an invertible change of variables because p is odd. Then Q = c² + uw. For each c:
- if 1 − c² ≢ 0 (p − 2 values of c), uw = 1 − c² has p − 1 solutions;
- if c = ±1 (2 values), uw = 0 has 2p − 1 solutions.

Total: (p−2)(p−1) + 2(2p−1) = **p² + p**.

**Lifting (a counted Hensel step).** Take a solution x mod p^j with j ≥ 1. Its p³ lifts mod p^(j+1) are u + p^j·t, where u is a fixed lift of x and t ranges over (Z/p)³. Since p^(2j) ≡ 0,

    Q(u + p^j t) = Q(u) + p^j · 2(u_c t_c + u_d t_d − u_b t_b)   (mod p^(j+1)),

so the lift is a solution exactly when t satisfies one affine-linear equation over F_p. Its coefficients (x_b, x_c, x_d) mod p are not all zero, because Q(x) ≡ 1. So exactly **p²** of the p³ lifts are solutions. Hence M = p^(2(k−1))·(p² + p) = p^(2k−1)(p + 1), and **a(p^k) = p^(2k−1)(p+1) + 2**. ∎

### p = 2
- **k ≥ 3:** this is A236554: a(2^k) = 2^(2k+2) + 32. The proof is in `../A236554/README.md`.
- **k = 1, 2:** direct counts give 8 and 64.

## Verification
- **Direct count.** `verify.py` counts the system directly for every n = 1..256, by grouping on a and taking exact convolutions of square counts. It compares the count with the closed form, then checks the closed form against all 49 OEIS data terms. Everything matches. The counter itself agrees with naive 4-nested-loop enumeration for n = 1..16.
- **Formal proof.** `lean/OeisLean/A236553.lean` (about 600 lines) proves:
  - `A_mul` (multiplicativity, via `ZMod.chineseRemainder`);
  - `A_odd_prime_pow` (the unit-or-zero split, the base count p² + p, and the counted Hensel step);
  - `A_two`, `A_four` and `A_two_pow` (the last reusing `A236554.card_A236554`);
  - `A_closed_form`, via `Nat.multiplicative_factorization`.

## Files
| File | What it is |
|---|---|
| `README.md` | This write-up: the result and the proof. |
| `layman.md` | Plain-English explanation. |
| `verify.py` | Direct count of A236553 and A227867 for n ≤ 256, compared with the closed forms and the OEIS data. |
| `bfile.py`, `b236553.txt` | b-file for n = 1..10000 from the closed form, checked against the direct count for n ≤ 256. |
| `lean/` | Lean 4 + Mathlib formal proof (`OeisLean/A236553.lean` together with `OeisLean/A236554.lean`, which it uses). Build with `cd lean && lake exe cache get && lake build`. |

## Log
- 2026-09-26: Derived the closed form (valuation split, p² + p points mod p, Hensel lifting). `verify.py` confirms it for n ≤ 256 and against all OEIS data.
- 2026-09-26: Formalized in Lean (`A_closed_form`); b-file written.
