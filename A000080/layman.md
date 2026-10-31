# A000080 in plain English

*Counting the leanest ways to connect points with triangles, and how a picture of the problem turned a
day-long computer search into a formula.*

---

## 1. The idea

Scatter some dots on a page. You may draw triangles, each one joining three of the dots. You want every
dot to be reachable from every other by hopping along triangles, and you want to be stingy: no triangle may be
spare. If you rub out any single triangle, something must fall apart.

Such a drawing is a **minimal triangle graph**. Two drawings count as the same if one is just a
relabelling of the other's dots. The question is how many different ones there are on n dots.

The ordinary-graph version, connecting dots with *lines* and no spare line, is the familiar question of
counting **trees**. Triangles behave differently. Two triangles may share two dots, and several triangles
can share one pair of dots. Rubbing one out still disconnects its "private" third dot, so each one is still
needed.

## 2. Where this shows up

* It is the simplest case of a *minimal connected hypergraph*, which shows up wherever connections come in
  groups (three-way interactions, triples of variables in a constraint system).
* The sequence is old. It was in Sloane's 1973 *Handbook of Integer Sequences* and was computed by
  R. Bowen in 1967.

## 3. What A000080 counts

On 3 dots there is one: a single triangle. On 4 dots there is again one: two triangles sharing two dots
(ABC and ABD). On 5 dots there are two: two triangles meeting at one dot (ABC, ADE), or three triangles on a
common pair (ABC, ABD, ABE). The counts go

    1, 1, 2, 4, 9, 19, 48, 117, 307, 821, 2277, 6437, 18634, ...

and the companion sequence A048781 adds up the triangles in all of them: 1, 2, 5, 13, 35, 86, ...

## 4. Why it looked hard

The numbers grow more than threefold per dot. The known method built every drawing and then threw away
duplicates with a graph-isomorphism program. Martin Fuller's well-tuned program in 2015 needed about a day to
reach 25 dots, where there are 1.5 billion drawings. Going further that way costs about 3.4× more per dot.

## 5. What we found

**A picture that explains every such drawing.** Draw one node for each dot and one for each triangle, and
join a triangle to its three dots. A drawing is minimal exactly when every triangle node, removed, splits this
network. That forces a simple shape. The drawing is a *tree* of pieces glued at dots, and each piece is either

* a **lone triangle**, or
* a **block**: some dots joined in a "no weak point" pattern, where several triangles may run along the same
  pair of dots, and every triangle pokes out one private corner.

Tree-like objects built from known pieces are exactly what a 1930s counting method (Pólya's) handles well.
With it the whole count becomes a formula, adding up over the ways to split a number into parts. No drawing
is ever built. All counts up to **80 dots** now take about two minutes:

* 26 dots: **4,987,668,376**
* 27 dots: **16,341,943,122**
* 80 dots: **2,617,646,724,739,221,521,059,589,519,172,655,683,909,550**

The same formula, with each triangle tagged, also extends A048781 (the total number of triangles).

You can browse examples at https://guillochon.github.io/minimal-triangle-graphs/: every drawing up to 9 dots,
and random ones beyond.

## 6. How we know it's right

* **It reproduces all 23 previously known counts** of both sequences, up to 25 dots. These came from a
  completely different method (search plus isomorphism tests).
* **The pieces match.** We took every drawing up to 15 dots, cut each one into its pieces, and counted the
  distinct pieces. They agree with the formula's own prediction of the pieces.
* **Two independent programs** (exact fractions in Python and fast modular arithmetic in Rust) give the
  same numbers.
* **Each big number was computed five times**, modulo five different large primes. Three determine it and
  the other two must agree, which catches arithmetic slips.
