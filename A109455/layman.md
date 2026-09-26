# A109455 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea: yes/no decisions by weighted voting

Picture a small committee deciding yes or no on a proposal. Each member votes yes or no. Some
members count for more than others, and the proposal passes if the weighted "yes" votes reach a
quota. Some examples:

- **Simple majority of 3:** each vote counts 1, and at least 2 are needed.
- **Unanimity:** each vote counts 1, and all of them are needed. (This is the logical AND.)
- **"Anyone can say yes":** each vote counts 1, and 1 is enough. (This is the logical OR.)
- **A chair with extra weight:** the chair counts 2, two ordinary members count 1 each, and 3 is
  needed. The proposal passes when the chair and at least one other member say yes. The two
  ordinary members together (1 + 1 = 2) can't pass it without the chair, so the chair has a veto.
  Different weights give genuinely different rules.

Mathematicians allow one more twist. A member's weight may be **negative**, so that their "yes"
counts *against* the proposal: "pass if Alice says yes and Bob says no."

Any yes/no rule that can be written this way (weights, a quota, pass if the weighted total reaches
the quota) is called a **threshold function**. Plenty of rules can't be written this way. The
classic example is "pass if exactly one of the two says yes" (XOR). No choice of weights and quota
produces it.

## 2. Where threshold functions show up

- **Voting systems.** Shareholder votes, the EU Council's qualified-majority rules, and the UN
  Security Council (whose veto powers can be expressed with weights) are all built this way. Game
  theorists call them *weighted voting games*.
- **Artificial neurons.** The simplest "neuron" in a neural network, the *perceptron*, takes
  yes/no inputs, multiplies each by a weight, adds them up, and fires if the total passes a
  threshold. That's exactly a threshold function. Asking "how many threshold functions are there?"
  is the same as asking "how many different behaviours can a single neuron have?". This question
  goes back to the 1960s, when it helped explain what one neuron can and cannot do.
- **Circuit design.** Engineers have studied building computers from threshold gates rather than
  ordinary AND/OR gates, and counting what those gates can compute was part of that work.

## 3. What A109455 counts

A109455(n) is the number of **genuinely different** threshold rules for a committee of n members,
where two rules count as the same if they differ only in **which member is which**.

For example, "pass if Alice says yes" and "pass if Bob says yes" are the same rule with the names
swapped, so they count once. With 2 members there are 14 threshold rules, but only **10** once name
swaps are ignored:

| | rule |
|---|---|
| 1 | always no |
| 2 | always yes |
| 3 | one particular member decides alone ("Alice decides" ≈ "Bob decides") |
| 4 | one member decides, but backwards (passes when they say no) |
| 5 | both must say yes (AND) |
| 6 | at least one says yes (OR) |
| 7 | not both say yes (NAND) |
| 8 | neither says yes (NOR) |
| 9 | one says yes and the other says no |
| 10 | one says yes, or the other says no |

So the sequence starts 2, 4, 10, 34, 178, 1720, … for committees of 0, 1, 2, 3, 4, 5, … members.

## 4. Why it's hard

The number of threshold rules explodes as the committee grows:

| members | all threshold rules (A000609) | different up to renaming (A109455) |
|---:|---:|---:|
| 6 | 15 million | 38,456 |
| 7 | 8.4 billion | 2,490,634 |
| 8 | 17.6 trillion | 550,112,272 |
| 9 | 144 quadrillion (1.4×10¹⁷) | 448,764,716,674 |

You can't list 144 quadrillion rules one by one and then group them. So the trick is to count the
groups **without listing their members**.

Our method rests on a classic counting idea (Burnside's lemma). Instead of counting groups
directly, you count, for every way of shuffling the members, how many rules *don't change* under
that shuffle, and then average. The rules left unchanged by a shuffle turn out to be much simpler
objects: rules on a small grid of points instead of on all possible voting patterns. Counting
those is far cheaper. The one really big term, "no shuffle at all", equals the total number of
threshold rules, which other researchers have already computed.

## 5. What we found

**The OEIS had a wrong number.** For 6 members the OEIS listed 590,440. That's impossible, and you
can see why without a computer:

- If you *also* allow flipping individual members' votes (treating "yes" as "no" for some
  members), there are only **1,113** different rules for 6 members (another OEIS entry, A000617).
- Taking away the flipping can split each of those 1,113 families into at most 2⁶ = 64 pieces.
- So there can be at most 1,113 × 64 = **71,232** rules, far fewer than 590,440.

The correct value is **38,456**.

**New terms.** Nobody had computed the values for 7 and 8 members. They are **2,490,634** and
**550,112,272**. For 9 members it is **448,764,716,674**. That run also recounted all 144 quadrillion rules
for 9 inputs from scratch, confirming a published number that had rested on a single 2006 thesis.
Along the way it found that a related OEIS entry (A000617) also seems to have a wrong value for 9
inputs, and it gave new values for two more entries (A002078 and A001529).

## 6. How we know the numbers are right

A new number in a reference work needs strong evidence, so we checked it several ways:

- **Two independent methods for 6 members.** Brute force (list all 15 million rules and group them)
  and the clever counting method give the same answer.
- **Rediscovering known results.** In the "no shuffle" case, the method recomputes numbers other
  researchers published long ago, exactly: all threshold rules for up to 8 members, plus two related
  counts. If the program had a bug, these would almost certainly come out wrong.
- **Divisibility check.** The averaging step divides a huge sum by n! (for example 8! = 40,320). A
  wrong sum would almost never divide evenly, and ours always does.
- **A certificate for every rule.** Each rule the program counts comes with explicit weights and a
  quota proving it really is a threshold rule. None failed.
- **Two implementations.** The Python version and the fast Rust version agree on everything they
  both computed.

## 7. Why it matters (honestly)

This is a small piece of mathematics, not a breakthrough. It still has real value:

- **Fixing an error in a widely used reference.** Researchers look up the OEIS to recognise
  sequences in their own work, and a wrong term can send someone down a dead end.
- **Extending what's known.** The new terms tell us how many essentially different weighted voting
  rules (or single-neuron behaviours) exist for committees of 7, 8 and 9 members.
- **The approach is what counts.** We didn't win with a bigger computer. We found a better way to
  count, and that is exactly what hard OEIS sequences usually need.
- **A clear frontier.** Going to 10 members would need the total number of threshold rules for 10
  inputs. That number is unknown, and computing it is a famous hard problem in its own right. The
  limit is now a clean, well-defined open question.

---

### Tiny glossary
- **OEIS:** the On-Line Encyclopedia of Integer Sequences, a huge reference catalogue of number sequences.
- **Threshold function:** a yes/no rule of the form "weighted sum of the inputs ≥ quota".
- **Up to renaming (permutation):** two rules are the same if relabelling the members turns one into the other.
- **Burnside's lemma:** a counting shortcut. The number of groups equals the average number of things left unchanged by each shuffle.
- **Perceptron:** the simplest artificial neuron. It computes exactly a threshold function.
