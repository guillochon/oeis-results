# A006455 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

Think of a to-do list with rules like "wash before drying" or "buy flour before baking". A set of
such "this before that" rules is called a **partial order**, or **poset**. Some pairs of tasks are
ordered and others aren't; you may fold laundry and bake in either order.

Now number the tasks 1, 2, …, n, and ask that the numbering never breaks a rule: whenever task x
must come before task y, x has the smaller number. A poset together with such a numbering is a
**naturally labeled poset**. Equivalently, it is a set of rules "i before j", always with i < j,
that is consistent: if 1 is before 2 and 2 is before 3, then 1 is before 3.

**A006455(n)** counts the naturally labeled posets on the tasks 1, …, n.

Tiny example, n = 3 (a(3) = **7**): no rules at all; one rule (1<2, 1<3 or 2<3); 1 before both
others; both others before 3; or the full chain 1<2<3. The pair "1<2 and 2<3" alone is not
allowed, because it forces 1<3 as well.

## 2. Where this shows up

- **Scheduling and sorting.** A numbering that respects the rules is a valid order to do the tasks.
  So the sequence also counts pairs (set of rules, valid schedule), up to renaming the tasks.
- **Matrices.** It also counts n×n upper-triangular 0/1 matrices with 1s on the diagonal that equal
  their own square in Boolean arithmetic (idempotent matrices).
- **Current research.** Recent papers study these objects via pattern-avoiding permutations
  (Bevan, Cheon and Kitaev, 2023), matrix methods (Cheon, Giraudo and coauthors, 2025–2026), and
  representation theory (Mazorchuk, 2026). Dividing by the count of all posets gives the average
  number of valid schedules of a random poset (Brightwell, Prömel and Steger, 1996).

## 3. What was known

The OEIS listed a(0) to a(12):
1, 1, 2, 7, 40, 357, 4824, 96428, 2800472, 116473461, 6855780268, 565505147444, 64824245807684.
The last three came from a program by Donald Knuth in 2001, and nothing new had appeared since.
A 2025 paper by Cheon, Giraudo, Kwon and Lee still says the sequence "is known up to n = 12".

## 4. Why it's hard

The numbers grow very fast; each term is about 100 times the one before. The natural method builds
the posets one task at a time: the new task n+1 goes last, and you choose which earlier tasks
must come before it. That visits every poset on the way, and a(13) alone is about 10 quadrillion
(10¹⁶).

## 5. What we found

**Four new terms:**

> a(13) = **10,252,978,223,661,714**
> a(14) = **2,223,739,823,604,074,838**
> a(15) = **657,800,677,264,030,983,838**
> a(16) = **264,134,770,940,774,037,129,920**

The key idea is to **cut the task list in two**: the tasks 1…n and the tasks n+1…n+k. Each half
is a naturally labeled poset on its own. What ties them together is, for each later task, the set
of earlier tasks that must precede it. Those sets have to fit together consistently, and counting
the consistent choices is a much smaller problem. So instead of building all 10¹⁶ posets for
a(13), we only need pairs of small posets: all 2045 shapes of 7-task posets combined with the 318
shapes of 6-task ones, each pair weighted by how many numberings it has. The computer then counts
the ways to glue each pair.

The gluing fact itself is a known theorem about posets (a Campo and Erné, 2018, building on
Erné's work from the 1970s on counting posets). What is new here is using it for naturally labeled
posets, plus a fast program for the gluing counts.

## 6. How we know it's right

- **Many cuts, one answer.** The same number can be computed by cutting at different places. For
  a(13) we cut 7+6, 6+7, 8+5 and 5+8; for a(14), 7+7, 8+6 and 6+8; for a(15), 8+7, 7+8 and 9+6; for a(16), 8+8, 9+7 and 7+9. The
  intermediate calculations are completely different each time, yet every cut gives exactly the
  same number, to the last digit. A bug would almost certainly break that agreement.
- **Old terms reproduced.** The same program recomputes every known term from a(2) to a(12),
  including Knuth's a(12), from ten different cuts.
- **The building blocks are checked.** The list of poset shapes comes from nauty's `genposetg`
  (Brinkmann and McKay's generator), and its counts match the known numbers of posets. From each
  shape the program reproduces two known sequences exactly for every n ≤ 9: A006455 itself and
  A001035 (all labeled posets).
- **The gluing formula is machine-checked.** The **Lean** proof assistant verifies every step of
  the cut-and-glue correspondence, and the resulting formula for a(n+k). It also verifies the fact
  behind the weights: how many numberings a poset shape has. A brute-force check of the formula on
  every cut up to 7 tasks is included too. What Lean does *not* check is the big computation
  itself; that rests on the agreement between different cuts.

## 7. Why it matters (honestly)

This is a solid but modest result, not a breakthrough.

- **A 25-year-old frontier moves.** The sequence is a classic (it appears in the original
  *Encyclopedia of Integer Sequences*), recent papers cite it, and it had been stuck at n = 12 since
  2001. Now it goes to n = 16.
- **Much less work.** Knuth's method visits every object. Cutting in half means the work depends
  on posets of about half the size, which is why a(15) took about a minute on a home PC and a(16) about 13 minutes.
- **It connects two traditions.** Erné's decomposition was used for counting *all* labeled posets
  (most recently to reach 19 points, Ayala 2026). It turns out to fit naturally labeled posets
  just as well.
- **Honest limits.** The key theorem is not ours, and the new terms rest on a computer
  calculation. It was checked many ways, but not proved by hand.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Poset (partially ordered set):** a set with consistent "this before that" rules, where some
  pairs may be unordered.
- **Naturally labeled:** numbered 1…n so that every rule goes from a smaller to a larger number.
- **Shape (unlabeled poset):** a poset with the names of the tasks forgotten; two posets have the
  same shape if renaming turns one into the other.
- **nauty:** a standard, widely used program package for listing combinatorial structures up to
  shape.
