# Square-filling motifs, in plain English

**Pictures:** https://guillochon.github.io/square-filling-motifs/

## The puzzle

Cut a square into an n×n grid of small tiles. Draw a path along the grid lines from the bottom-left corner to the
bottom-right corner, with one rule: every tile must "own" exactly one of the segments of the path, a segment on
one of its four sides. Since there are n² tiles, the path has exactly n² steps. Doug McKenna calls such a path a
**motif**.

Now replace every tile by a shrunken copy of the whole motif, turned so that it runs along that tile's own segment
in the same direction. The result is a motif on an n²×n² grid. Do it again and again, and the path wiggles ever
more finely until, in the limit, it passes through every point of the square: a space-filling curve, like
Hilbert's famous one.

The catch: the copies sit side by side, and their paths can bump into each other. For almost every motif, after one
or two rounds two pieces of path land on the same grid line or meet at a point. McKenna calls the rare motifs that
never bump into themselves, however many rounds you do, **globally self-avoiding**. OEIS
[A350148](https://oeis.org/A350148) counts them, one per mirror pair, on the 1×1, 3×3, 5×5, … grids.

## What was known

On the 5×5 grid there is exactly one (McKenna found it in 1979). On 7×7 there are 7. On 9×9 there are 10,101, and on
11×11 there are 20,305,328; that last count took McKenna 147 hours of computer time, listing every motif one by one.
For 13×13 he found examples but no count.

## What is new

| Grid | Motifs |
|---|---|
| 13×13 | 7,287,109,207,271 |
| 15×15 | 3,892,817,988,876,823,083 |
| 17×17 | 370,720,354,045,169,062,282,196,655 |
| 19×19 | being computed |

## How

McKenna proved two remarkable facts about these motifs. First, if you glue the square's opposite edges together
(making a doughnut), a good motif becomes a single loop that passes through every grid point exactly once. Second,
the ring of tiles around the edge of the square is completely determined by a short code: a string of 0s and 1s,
the same read forwards or backwards, with no two 1s side by side. For a 13×13 grid there are only five such codes.

With the outer ring fixed, the rest is a self-contained puzzle: give every inner tile one of its sides so that the
chosen sides form a single path through every inner grid point. Instead of listing the answers, a computer can
count them by sweeping across the grid one tile at a time, remembering only how the path crosses the current
"front" line (which loose ends will later join up). Many different partial drawings look the same at the front, so
they are counted together. This counts all 20 million 11×11 motifs in a tenth of a second, 13×13 in eleven
seconds, and the 17×17 count, about 3.7·10²⁶, in an hour and a half.

## How we know it is right

- The method reproduces every count McKenna published, including how his 9×9 and 11×11 motifs split between the
  edge codes.
- An independent test, written straight from the definition, iterates each motif and checks that no copies ever
  touch. Every motif up to 9×9 passes it, as do random samples of the larger ones.
- A short argument plus a computer check shows that every motif the method counts must pass that test. Copies of a
  motif only meet along their outer rings, and those are fixed by the code, so it is enough to check how the rings
  meet in every way they can.
