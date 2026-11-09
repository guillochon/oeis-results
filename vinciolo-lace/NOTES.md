# Vinciolo lace: working log

The self-contained write-up is `PUBLIC.md` (published as README.md in oeis-results/vinciolo-lace). This file
is the log of how we got there, published as NOTES.md.

## 2026-11-07

* Started from a photo of page 37 recto of Vinciolo's *Les Secondes Oeuvres* (1603), "Carrez contenant 41 maille",
  posted online. Question: is it related to the symmetric sandpiles of A348567? Traced the plate (blob lattice
  walk, 107 × 82 cells); both central panels are 43 × 43 and nearly fully symmetric (lower 94%, upper 85%).
  Compared with the 43 × 43 sandpile identity under seven height thresholds: best agreement 55–59% against 50% for
  random symmetric patterns. Not related; the panels are not recurrent configurations either.
* Tried rules that both panels satisfy:
  * open meshes 8-connected (both panels; random symmetric designs only 1–5% of the time). Symmetric count
    1, 6, 43, 556, 14125, 721874, 74784703 (odd n ≤ 13; even n equals n − 1); without symmetry this is A286139.
    Random designs under it look like noise.
  * no lone cell (8 neighbors of the other color): the lower panel breaks it (flower centers).
  * outward growth from the center, single seed or per piece: the panels break it about as often as random designs.
  * "only 3 × 3 (3 × 4, 4 × 4) blocks Vinciolo uses": zero violations by construction and lace-like samples
    for 3 × 4, but this was fitted to two panels only.
* Downloaded the Met's 1603 pages (51 web-size images, then 28 full-resolution originals, 72 MB; public domain).
  Traced all; 14 of the 41-mesh plates without lattice conflicts. Panel finder: square windows that agree with
  their D4 images, ranked by size × symmetry, required to sit on the solid frame (outer ring ≥ 90% printed);
  without the frame condition two windows cut through designs and created false pinholes at their edges.
  19 panels (13 pages), symmetrized by majority vote. The Met trace of page 37 agrees with the photo trace on 98%
  of cells (99.6% for the symmetrized upper panel).
* With 19 panels the two-panel rules all fail (every new panel uses unseen 3 × 3 blocks); all 512 blocks occur
  except the pinhole (open cell, all 8 neighbors printed, frame printed): expected about 120, observed 0. Next
  rarest block: 16 uses in one panel. 3 × 4 and 4 × 4 seen-only rules fail leave-one-out. Open-region
  connectivity holds in 10/19; lone printed dots in 7/19 panels.
* Sequence: D4-symmetric n × n, no pinhole. Color swap (1 = open) gives "every 1 has a king neighbor 1, outside 0",
  the symmetric analog of A202788 (brute force reproduces A202788 for n ≤ 5). Brute force n ≤ 14.
* Transfer matrix over the octant (`tools/pinhole_tm.cpp`), 2 bits per frontier cell, h + 1 frontier cells.
  Exact 384-bit counts to n = 43 (n = 45 ran out of WSL's 23 GB), then residues mod 2^64 and four primes near
  2^62 for n = 43..47 (`crt.py`), exact and CRT agreeing at n = 43. a(1)–a(47) in `b_pinhole_d4.txt`.
* PARI brute force (`tools/check.gp`) agrees for n ≤ 10. Draft submission `oeis_submission.md`.

## 2026-11-08: growth and region rules

Asked whether the panels follow any rule beyond "no pinhole", in particular one about how the design grows.
Tools: `tools/growth_rules.py` (rule battery), `tools/small_open.py`, `tools/null_sampler.c` (20 random
no-pinhole designs per panel, same size and darned fraction, `null/`), `tools/count_minregion.c`.

- **Radial growth rules all fail, in every panel.** Tested: open or darned cells must have a same-color neighbor
  closer to the center (or, the other way, closer to the edge), except the seed ring of their own region. Three
  radii (square, diamond, round) and 4- or 8-neighbors, 24 variants. Every variant has violations in all 19 panels.
  The motifs bend back toward the center (scrolls, hooks, flower heads), so neither color grows monotonically.
- **No open region has fewer than 5 meshes.** This holds in 19/19, and size 5 occurs once: the open cross at the
  center of p23r. Regions are 8-connected open meshes, and the outside counts as darned. No pinhole is the case
  "size at least 2". Random designs of the same size and density break the rule in most panels. They pass at a rate
  of 0% on p23r and p24r (the dense ones) and 25–90% elsewhere, so the product over the 19 panels is below 1e-8.
  The reason is that Vinciolo's open meshes form one big web: 1 region in 10 panels and at most 10 in the others.
  Random designs have a median of 2–91 regions, most of them tiny.
- Darned islands (not touching the frame) of size 2 occur in p26r_0 and p33r_1, so there is no matching rule on
  the darned side. Lone darned dots (flower eyes) are common.
- Counts of the generalized rule, "every 8-connected region of 1s has at least k cells", for k = 1..8, are in
  `tools/minregion_counts.txt`. k = 2 reproduces a(1)–a(12). For k = 5:
  - symmetric: 1, 1, 5, 6, 42, 43, 601, 612, 17134, 17756, 990833, 1036978;
  - all arrays: 1, 1, 227, 36632, 15691156.
  Neither the k = 3, 4, 5 families nor their unsymmetric versions are in the OEIS. The unsymmetric k = 2 version
  is A202788. Going further needs a transfer matrix that tracks region labels and sizes up to k − 1.

## 2026-11-08: the minimum-region table T(n,k)

T(n,k) = number of D4-symmetric n x n 0/1 arrays in which every 8-connected region of 1s has at least k cells
(k = 2 is the no-pinhole sequence, and k = 5 is the strongest region-size rule all 19 panels obey). The user chose
a single table entry read by antidiagonals. Draft: `oeis_submission_table.md`; b-file: `b_minregion_table.txt`
(antidiagonals n + k = 2..34, 561 terms).

- `tools/minregion_tm2.cpp`: a transfer matrix over the octant. A quotient component touches the diagonal wall
  and/or the vertical-midline wall when one of its cells is equal or king-adjacent to its own mirror image. Its
  full-array region then has stabilizer <those reflections> (order 1, 2, 2 or 8), and size
  (sum of orbit sizes) * |H| / 8. States are packed into 192 bits; counts are kept mod 2^128 and mod 2^61 - 1
  and combined by CRT (exact up to 2^189). `minregion_tm.cpp` is the first, string-keyed version.
- Checks: brute force (`count_minregion.c`) n <= 12, k <= 8, 96 values; the two TM versions agree; column 2
  matches `b_pinhole_d4.txt` through n = 32 (past 2^128, so the CRT path is exercised); PARI
  (`tools/check_table.gp`) agrees on antidiagonals 2..10; `assemble_table.py` checks monotonicity in k and
  T(n,k) = 1 for k > n^2.
- Runs: `tools/run_table.py` (resumable; `table_runs.txt`, logs `table_log1..3.txt`). 4 workers x 4 GB through
  d = 32, then 2 x 9 GB, then 1 x 14 GB for (29,5) and (27,7) at about 50 million states each. Antidiagonal 35
  needs bigger jobs; stopped.
- Not in the OEIS (Nov 8): the table, columns 3, 4, 5, 6, 10 and rows 5 and 6. k = 1 is 2^A008805(n-1).

## 2026-11-08: patch-mix (maximum-entropy) designs

The user asked for a page section that uses the distribution of 3 x 3 blocks Vinciolo used.

- `tools/blockstats.py` -> `block_target.txt`: counts of the 102 symmetry classes over the 47,867 windows (frame
  darned, as in build.py); 101 are used, and the unused one is the pinhole.
- `tools/block_fit.c`: P(design) ~ exp(-sum of w[class(window)]) over D4-symmetric framed panels, with the pinhole
  forbidden. Fitted by persistent-chain moment matching, one chain per panel size (`block_sizes.txt`), 2500
  iterations, with the weights averaged over the last quarter. Output: `block_weights.txt`, `block_fit_classes.txt`,
  `block_samples_43.txt`, log in `tools/block_fit_log.txt`.
- Fit, from fresh chains: KL(Vinciolo || model) = 0.040. Density-matched rule-only designs give 0.207. Patch kinds
  within ±50% of Vinciolo's rate cover 87% of his windows for the model and 59% for rule-only designs (from
  `tools/block_compare.py` -> `block_compare.json`). The model stalls near KL 0.03-0.04: the chains switch between
  dense and open phases during the fit, and neither a single size nor per-panel sizes did better.
- The samples look like lace texture (bars, frames, stepped diagonals, medallions), but they have no large motifs.
  The page's JS generator is checked in Node: KL 0.042 against the data, no pinholes.

## 2026-11-09: off-center pocket rule

- The only open region under 9 meshes in the 19 panels is the 5-mesh cross at the center of page 23. Leaving out the
  region through the center cell (2 x 2 block for even n), the smallest region is 9: eight mirror copies in page 31
  (next: 66, four copies in page 24; then 192 and up). 25 off-center regions in all. `tools/offcenter_pockets.py`.
- Density-matched no-pinhole nulls (null/): joint pass chance for "off-center min >= 9" ~1e-11, for ">= 5" ~3e-9
  (zero-pass panels counted as 0.5/20, so rough). Rests on one panel for the 9.
- Not added to the T(n,k) table (decided with the user). Used instead as the rule for the site's random designs:
  editor default k = "9, except through the center"; patch-mix generator rejects flips that break it, starting from
  all open. `tools/pocket_patchmix.js` runs the page's own generator in Node: with the rule KL 0.050, 86% within
  +-50%, 57/57 designs obey it (independent flood fill), 115 ms per design; without the rule KL 0.044, 25/57 obey.
  Writes block_fit_classes_rule.txt (now used by block_compare.py) and block_samples_43_rule.txt (page gallery).
