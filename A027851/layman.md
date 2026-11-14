# A027851 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

A multiplication table tells you how to combine any two things from a list into a third. Ordinary
multiplication has one property we rarely think about: when you multiply three numbers it does not matter
which pair you do first, (a × b) × c = a × (b × c). A table with that property is called a **semigroup**.
It is about the least you can ask of a "multiplication": no identity, no division, not even a × b = b × a.

Gluing words together is a semigroup ("ab" then "c" is the same as "a" then "bc"). So is "do this step,
then that one", which is why semigroups describe the behavior of machines and programs.

## 2. Where this shows up

- Every finite-state machine has a semigroup that records what sequences of inputs do to it.
- Groups (symmetries) and the sets of matrices used throughout physics and computing are semigroups with
  extra structure.
- How many semigroups there are of each size is one of the oldest counting questions in algebra, started
  by hand in the 1950s and then one of the first to be handed to computers.

## 3. What A027851 counts

The number of different semigroups with n elements, where two tables are the same if one becomes the other
by renaming the elements.

With 2 elements there are 5. The sequence begins 1, 1, 5, 24, 188, 1915, 28634, 1627672, …

There is a second, coarser way to count: also call two tables the same if one is the *mirror image* of the
other (read every product backwards, b × a in place of a × b). That count is the sequence A001423.

## 4. Why it's hard

The numbers explode: with 10 elements there are more than 10^19 semigroups. Nobody lists them. In 2012 a
team in St Andrews found the mirror-image count for 10 elements, 12,418,001,077,381,302,684, with a clever
formula for the overwhelming majority and a large distributed computer search for the rest.

The plain count (no mirror images) was left open. To get it from the other one you need to know how many
semigroups **are their own mirror image**: each of those was counted once, and every other semigroup was
counted once for a pair. The self-mirror ones are rare, one in twenty million, and hidden among the rest.

## 5. What we found

We counted the self-mirror semigroups with 10 elements directly:

**613,649,634,393**

and therefore the number of semigroups with 10 elements is

**24,836,001,541,112,970,975.**

The key idea: a semigroup is its own mirror image exactly when some renaming of its elements turns the
table into its mirror. So instead of searching all tables and testing each, we fix the renaming first and
build only tables that respect it. Half of every table is then forced by the other half, and semigroups
without the symmetry are never looked at. A standard counting rule (Burnside's lemma) turns the totals for
the 42 kinds of renaming into the answer. It took 35 minutes on a desktop PC.

## 6. How we know it's right

- **Smaller sizes.** The same program gives the known answers for every size up to 9, including the
  breakdown by number of idempotents (elements with a × a = a). That tests the whole method.
- **A second road to the same partial counts.** Semigroups in which a × b = b × a are their own mirror
  image, and their number for 10 elements was published by another author with another method. Our program
  finds exactly that number among its results. The same holds for a second family (inverse semigroups).
- **Whole numbers.** The method adds up fractions with denominators as large as 3,628,800; every total is
  a whole number.
- **A caveat.** The final count uses the 2012 mirror-image count. If that were wrong ours would be too; the
  self-mirror count would not.

## 7. Why it matters (honestly)

This fills a known gap rather than opening new ground.
- It completes a classical reference sequence one step further, a step that had been open since 2012.
- It shows the gap needed an idea rather than more computing power: 35 minutes, against a
  distributed computation for the 2012 count.
- Ten related entries (semigroups by number of idempotents, bands) get their order-10 values at the same time.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number
  sequences.
- **Semigroup:** a multiplication table where (a × b) × c = a × (b × c).
- **Self-dual (self-converse):** a semigroup that is the same as its mirror image, up to renaming.
- **Idempotent:** an element with a × a = a.
- **Burnside's lemma:** a rule for counting things up to symmetry by averaging over the symmetries.
