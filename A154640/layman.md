# A154639–A154645 in plain English

*Packing ever-smaller balls into the gaps between balls, and counting how many arrive at each step.*

---

## 1. The packing

Put three marbles of the same size flat on a table, all touching each other, and lay a sheet of glass on top so it
touches all three. The two flat surfaces and the three marbles are five objects that all touch one another.

Now fill the gaps. Any four of the five objects leave room for exactly one more ball that touches all four.
Add it, and repeat with every new group of five touching objects. The balls get smaller and smaller and fill the
space between the plates, like a three-dimensional version of the famous Apollonian circle picture.

## 2. Generations

Each round of filling is a **generation**: generation 0 is the starting five, generation 1 adds 5 balls,
generation 2 adds 20. In 2009 Colin Mallows counted the balls in each generation up to the sixth, along with
the sum of their curvatures (curvature = 1/radius).

There was a snag. In generation 2, some of the new balls touch each other, forming new groups of five touching
balls earlier than expected. Should balls built from those groups count in generation 3 or generation 4?
Mallows counted both ways, calling them strategy (a) and strategy (b). That gave six OEIS sequences, plus a
seventh counting the "words" that describe the steps.

## 3. What is new

- **Exact formulas** for strategy (a): the number of balls in every generation, and the sums of their curvatures
  and squared curvatures, are given by simple fractions of polynomials. Behind this is a theorem about
  reflection groups: each ball corresponds to exactly one coset of a smaller group of reflections.
- **Mistakes fixed.** Six of Mallows's numbers were wrong. The clearest is generation 6 in strategy (a): he wrote
  3330 balls, but his own curvature totals for that generation fit exactly 2340, the correct count.
- **Ten more generations** for strategy (b), up to generation 16, where more than six billion balls are added.
  Every ball is tracked with whole numbers, so there is no rounding error.
- **A guess** that strategy (b) also follows a simple formula, true for all 17 known generations.

## 4. See it

The picture page shows the packing in 3D, colors each ball by generation, and highlights the balls that
strategy (b) counts earlier: [Apollonian Sphere Generations](https://guillochon.github.io/apollonian-spheres/).
