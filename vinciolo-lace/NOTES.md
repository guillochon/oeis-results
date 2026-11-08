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
