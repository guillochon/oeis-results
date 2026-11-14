# A132183 in plain English

*What this sequence counts, why anyone cares, and what we found.*

---

## 1. The idea

A committee of n people votes yes or no. A **voting rule** says which groups of yes-voters are enough to
pass a motion. Two things should be true of any sensible rule:

* **More support never hurts.** If a group can pass a motion, so can that group plus anyone else.
* **Some members matter more.** Suppose the members can be ranked, from the most influential (the chair)
  down to the least. Then swapping a member of a winning group for a more influential one never turns a
  win into a loss.

A rule with both properties is called a *complete simple game* by people who study voting, and a *regular
Boolean function* by people who study circuits. "The motion passes if the yes-voters hold more than half of
the shares" is one example, for any way of handing out shares. But there are many regular rules that no
share count can describe.

## 2. Where this shows up

* **Voting power.** Councils, boards and parliaments with weighted votes are regular rules. To ask "which
  rule gives each member a fair amount of power?" one has to search the list of all rules.
* **Circuits and threshold logic.** A neuron-like gate that fires when a weighted sum passes a threshold is
  a regular function. Regular functions are the larger, easier-to-describe family that contains all of them.
* **Knuth's The Art of Computer Programming.** Volume 4A studies the pattern behind these rules (the
  "binary majorization lattice") and gave the counts up to 9 variables in 2007.

## 3. What A132183 counts

The number of regular rules for n ranked members, including the two trivial rules "always pass" and "never
pass".

For 2 members (A stronger than B) there are 5: never pass; pass only if both say yes; pass if A says yes;
pass if anyone says yes; always pass. Note that "pass if B says yes, whatever A says" is not on the list:
A is the stronger member, so A alone must be at least as good as B alone.

The sequence starts 2, 3, 5, 10, 27, 119, 1173, 44315, 16175190, 284432730176.

## 4. Why it's hard

The counts explode. Going from 8 to 9 members multiplies the count by 17,584; from 9 to 10, by about 15
million. In 2010 a research group listed all 284 billion rules for 9 members one by one, which took
44 minutes. Listing the rules for 10 members the same way would take about a thousand years, and for 11
members far longer than the age of the universe. So the next terms needed a way to count the rules without
looking at each one.

## 5. What we found

Two new terms:

* **a(10) = 4,191,080,125,665,004,647** (about 4 quintillion)
* **a(11) = 7,950,574,918,162,482,690,174,494,937,400** (about 8 nonillion)

The key idea: split the committee into its 4 strongest members and everyone else. A rule for the whole
committee is then 16 smaller rules for "everyone else", one for each way the top 4 can vote, and these 16
rules must nest inside each other in a fixed pattern. Counting nested families of small rules can be done
with running totals instead of lists. For 11 members the small rules are the 44,315 rules for 7 members, and
the whole count takes a few minutes on a desktop PC.

## 6. How we know it's right

* **The old terms come out right.** The same programs give every known value up to 9 members.
* **Two different splits agree.** For 10 members we split off the top 2 (working with the 16 million rules
  for 8 members) and, separately, the top 4 (working with the 1,173 rules for 6 members). Both give the same
  19-digit number.
* **A side count matches other people's work.** The same tables also count the "self-dual" rules, where a
  group wins exactly when its complement loses. Those numbers were computed by others up to 11 members with
  a completely different technique, and ours match all of them, including the 16-digit value for 11.
* **A second program.** The 11-member count was recomputed by a separately written program that organizes
  the final sum differently. It agrees.

## 7. What's still open

The count for 12 members (a number with about 50 digits) is out of reach of this method on one PC. It needs
another idea.
