# A000112 in plain English

*Counting the shapes of "this before that" rule sets on 17 items, by borrowing a count of the labeled
ones.*

---

## 1. The idea

Think of a to-do list with rules like "wash before drying" or "buy flour before baking". A consistent set of
such "this before that" rules is called a **partial order**, or **poset**. Some pairs of tasks are ordered
and others aren't: you may fold laundry and bake in either order.

Now forget what the tasks are and keep only the pattern of rules. Two rule sets have the same **shape** if
renaming the tasks turns one into the other. **A000112(n)** counts the shapes on n tasks.

Tiny example, n = 3 (a(3) = **5**): no rules at all; a single rule with one task left over; one task before
both others (a "V", drawing earlier tasks lower); both before one (a "Λ"); and the full chain 1 < 2 < 3.

A companion sequence, A000608, counts the shapes that don't fall apart into separate pieces. It follows
from A000112 by a standard formula.

## 2. Where this shows up

- **Order and scheduling.** Posets describe dependencies in projects, builds and computer programs, and the
  "happened before" relation in distributed systems.
- **Mathematics.** Posets are the basic objects of order theory and lattice theory. On finite sets they are
  the same thing as certain topologies.
- **A classic.** The sequence was in Sloane's 1973 *Handbook of Integer Sequences*. Counting posets is a
  long-running benchmark for computer enumeration.

## 3. What was known

The OEIS listed a(0) to a(16):

    1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284, 46749427, 1104891746, ...

and a(16) = 4,483,130,665,195,087. a(15) and a(16) came from Brinkmann and McKay in 2002, who generated
every shape one at a time. Nothing new had appeared since.

Meanwhile the **labeled** count, where the tasks keep their names (A001035), kept moving. It reached 19
tasks in 2026 (Ayala). Labeled counts can use clever cut-and-glue formulas; shapes are harder, because you
must avoid counting the same shape twice under different names.

## 4. Why it's hard

The numbers grow fast, about 80 times per extra task near n = 17. Listing every shape and checking for
duplicates, as in 2002, would mean building **378 quadrillion** shapes. That is far beyond a home computer.

## 5. What we found

> a(17) = **378,180,373,015,320,248**
>
> and for the connected ones, A000608(17) = **373,626,005,846,478,701**.

The key idea is an old counting trick, **Burnside's lemma**. It says: the number of shapes equals the
average, over all ways of shuffling the task names, of the number of labeled rule sets that the shuffle
leaves unchanged.

For n = 3: there are 6 shuffles. Doing nothing leaves all 19 labeled rule sets alone. Each of the 3 swaps
of two tasks leaves 3 alone, and each of the 2 rotations of all three leaves 1 alone. The average is
(19 + 9 + 2) / 6 = **5**.

The useful surprise is what a shuffle can leave unchanged. Tasks that the shuffle cycles around can never
be ordered among themselves, and the cycles behave like single "super tasks". When the cycle lengths share
no common factor, as in the swaps and rotations above, the count is simply a **known labeled count on fewer
tasks**. For n = 17 those easy shuffles make up **99.34%** of the total, and they come straight from
A001035.

The remaining 0.66% comes from shuffles with several cycles of related lengths (mostly several swaps at
once). For those we computed a handful of exact counts on 11 to 15 labeled tasks: how many rule sets have a
particular small pattern in them. Some came from cutting each rule set into three blocks and gluing (the
same trick that extended A006455). Others came from visiting all 34 billion shapes on 13 tasks, which is
small enough to do.

The whole computation took a few days on a home PC. The final step that combines everything now runs in
about 1.6 hours.

You can see the method worked out on small cases, and every poset up to 7 tasks, at
https://guillochon.github.io/seventeen-point-posets/.

## 6. How we know it's right

- **Old terms reproduced.** The same program recomputes every known term a(1) to a(16). a(15) and a(16) go
  through the same hard-case machinery as a(17).
- **It divides evenly.** Burnside's sum must be a multiple of 17! (about 356 trillion). A wrong piece would
  almost never give a whole-number answer, and the program checks it.
- **Hard pieces computed twice.** Every exact count that could be computed by both methods (cut-and-glue,
  and visiting all shapes) was, and the two agree to the last digit.
- **Small cases by brute force.** The cycle trick was checked by brute force on up to 6 tasks.
- **The companion checks out.** The formula linking A000112 and A000608 reproduces all known terms of
  A000608, and runs backwards to give A000112 again.

## 7. Why it matters (honestly)

- **A 24-year-old frontier moves by one.** a(16) dates from 2002. Now there is a(17).
- **A cheaper route.** Instead of building all 378 quadrillion shapes, the method leans on labeled counts
  that other people had already computed, plus a few exact pattern counts.
- **Honest limits.** Burnside's lemma and the labeled counts are not ours. The new term rests on a computer
  calculation, checked many ways but not proved by hand. a(18) would need pattern counts on 16 tasks, a job
  on the scale of Ayala's labeled computation rather than a home PC.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalog of number sequences.
- **Poset (partially ordered set):** a set with consistent "this before that" rules, where some pairs may
  be unordered.
- **Labeled vs. shape:** a labeled poset keeps the task names; a shape (unlabeled poset) forgets them.
- **Burnside's lemma:** a classical rule that counts objects up to symmetry by averaging what each symmetry
  leaves unchanged.
- **nauty / genposetg:** a standard program package that lists posets one per shape (Brinkmann and McKay).
