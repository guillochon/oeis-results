# A350148: Hilbert-style square-filling curve motifs (McKenna)

Pictures: **https://guillochon.github.io/square-filling-motifs/** · plain English: [layman.md](layman.md)

a(n) counts the distinct globally self-avoiding motifs on the (2n+1)×(2n+1) subdivision: threads of tiles
from one bottom corner to the other whose iterated edge-replacement compositions G<sup>k</sup> never self-contact.
Known: 1, 0, 1, 7, 10101, 20305328. McKenna (ECA 2:2 (2022) #S2R13) enumerated a(5) by backtracking in 147
CPU-hours and, for 13×13, found only sample motifs from partial runs.

## Result

| n | grid | a(n) | per zipping mode | time |
|---|---|---|---|---|
| 2 | 5 | 1 | 1 | – |
| 3 | 7 | 7 | 7 | – |
| 4 | 9 | 10101 | 660 + 9441 | – |
| 5 | 11 | 20305328 | 58936 + 19854452 + 391940 | 0.1 s |
| 6 | 13 | **7287109207271** | 34174087 + 344435015915 + 3695591785 + 18614455315 + 6920329970169 | 11 s |
| 7 | 15 | **3892817988876823083** | 30561568363 + 5580524782288346 + 636529366875968 + 11306407057702863 + 3780168885982901976 + 10659326850390 + 63424610454389340 + 31690341344245837 | 9 min (peak 2.7M states) |
| 8 | 17 | **370720354045169062282196655** | 135879677410959 + 411297300212654023648 + 51711377581144855936 + 83957935550425652760350 + 16463131878906415553612712 + 41124585630148011913 + 211265298320859583865522 + 192143246383015915029769 + 796636287707588212 + 348048211035339330806375 + 306467364859207877763 + 856636468990701460108647 + 352564359608581854246244849 | 1.5 h on 6 threads (`big`; peak 40M states) |
| 9 | 19 | running (`bash run_n.sh 19 big 5 6`, results/n19.txt), 21 modes, roughly 2 to 3 days | | |

Modes are listed in McKenna's order (mode bits 0…01 first). The 9×9 and 11×11 per-mode splits are McKenna's,
reproduced exactly.

## Method

Two theorems from McKenna's paper turn the global condition into a local one:

- **Lemma 3.17.** A motif's tile path is a Hamiltonian cycle on the n×n torus. In the plane: every interior
  lattice point is visited, and of each pair of identified boundary points at most one.
- **Theorem 3.1 (Fibbinary zippers).** All 4(n−1) perimeter tiles are frozen by one of F<sub>(n−3)/2</sub> zipping
  modes, each a palindromic Fibbinary bit string.

`motif.boundary` turns a mode into the frozen perimeter. Each side, walked in the path's flow direction, starts on
the outer track and switches track at every 1 bit. It reproduces the E motif's boundary, and it must be the frozen
boundary for every mode listed, because the DP finds globally self-avoiding motifs with it (Theorem 3.1 makes the
boundary a function of the mode).

With the perimeter fixed, `rust/` counts by a cell-by-cell **frontier transfer matrix**. Every free tile owns one
of its four edges. No edge can be owned twice. The owned edges must form one path a → b through every interior
lattice point. The state holds the frontier points' degrees, strand labels and one "owned" bit per frontier edge.
11×11 peaks at 14k states and 13×13 at 190k.

For 17×17 and up, `rust/src/bin/big.rs` runs the same transitions (it calls the same `step`) with:
- 64-bit keys. Strand ends are written as brackets, two bits per frontier point, which works because the
  pairings below a frontier never cross.
- 128-bit counts.
- The next layer split into 256 hash shards, built on several threads and, when it would not fit the memory
  budget, a group of shards at a time, with the current layer spilled to disk.

It reproduces every per-mode count up to 15×15 with the same peak state counts, and the 17×17 mode 000000001
computed by the old engine (121 s against 778 s). The dominant 17×17 mode (001010101, 95% of a(8)) was rerun on the old engine and agrees
(`results/n17_check_old_engine.txt`, 2946 s).

## Checks

- **Known terms.** The DP reproduces every known term, mode by mode.
- **Exact test** (`motif.globally_self_avoiding`). This decides "G<sup>k</sup> self-avoiding for all k" without the
  paper's lemmas: it closes the set of 2×2 supertile patterns that occur in any G<sup>k</sup> under substitution and
  checks each block. It accepts McKenna's E motif and rejects his J. Every DP motif for 5×5, 7×7 and 9×9 (mode
  00001) passes, as do 40 random motifs per mode for 11×11 and 30 per mode for 13×13 (`verify.py`, uniform
  sampling from the DP).
- **Soundness for all counted motifs** (`sound.py`). Copies of G in neighbouring supertiles touch only through
  perimeter tiles. So if every locally plausible 2×2 tile pattern is safe for the frozen boundary, every DP motif
  is globally self-avoiding, by induction on k. A plausible pattern has distinct edges, at most one arrival and
  one departure per point, and parity changing exactly across owned edges (McKenna Lemma 3.6). All 192 plausible
  patterns are safe for every mode up to 19×19 (`results/sound_n19.txt` for the 21 modes of 19×19). The same test flags 164 and 186 of them for the Hilbert curve and
  McKenna's J.

Completeness rests on McKenna's Lemma 3.17 and Theorem 3.1. Soundness is the computed check above.

## Files

- `motif.py`: conventions, composition, exact global test, modes, frozen boundary
- `make_bounds.py`: writes `bounds/n{N}_{mode}.txt`
- `rust/`: `count FILE`, `enum FILE`, `sample FILE K SEED`
- `verify.py`: exact-tests DP output; `sound.py`: the boundary soundness check
- `results/`: run logs
