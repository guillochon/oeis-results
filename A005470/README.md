# A005470 — unlabeled planar graphs (and A003094, A021103)

**Status:** done to n = 30 for all three sequences (`tower_30.json`, b-files `b005470.txt`, `b003094.txt`,
`b021103.txt`, n = 0..30), checked as described under Checks.
**Previously known:** A005470 and A003094 a(0)..a(15) (a(13)-a(15) by Jeff Sponaugle, Sep 2026), A021103
a(0)..a(14), by exhaustive generation (plantri / nauty), and A021103(12)-(14) also by the formula of Gagarin,
Labelle, Leroux and Walsh (2009).
**Pictures:** the Planar Graph Atlas, https://guillochon.github.io/planar-graphs/.

## Results (n = 16..30; the full tables are the b-files)

| n | A005470 (all) | A003094 (connected) | A021103 (2-connected) |
|---|---|---|---|
| 15 | (known) | (known) | 1115272493022 |
| 16 | 52408033776878 | 49760001456551 | 22493653916620 |
| 17 | 1114782893212983 | 1059703295702139 | 462624367146705 |
| 18 | 24140058638845459 | 22969738086867520 | 9676037828904898 |
| 19 | 530861925075735556 | 505542319715022144 | 205363624344678050 |
| 20 | 11832495146001556996 | 11276122690200599855 | 4415099749028218062 |
| 21 | 266892197689888343244 | 254499290386523690634 | 96009973957009506442 |
| 22 | 6084033586294190154202 | 5804661393043936816211 | 2109220961695218583498 |
| 23 | 140012314430574390754809 | 133647000805745363204565 | 46763841586042636362681 |
| 24 | 3249784959167329158913333 | 3103364947739587802644026 | 1045438229538044361911334 |
| 25 | 76016793864553053210040815 | 72619635686723255869114898 | 23548181230365705962266586 |
| 26 | 1790733466003422184670125640 | 1711297840974659008783320156 | 534073930971881285532076914 |
| 27 | 42457813673889602379288800530 | 40587146973195981103012285191 | 12189355449099428544303427991 |
| 28 | 1012653135463833979592320804807 | 968313129216044117167721499740 | 279818237785356119481685635143 |
| 29 | 24285021128825709709135132968024 | 23227758888209768266519857560467 | 6457929182734030291113816963729 |
| 30 | 585345393678411726175403670215562 | 559996781949783047710910447811156 | 149782224970101334436971067029495 |

## How to reproduce
1. `bash run_gpu.sh 30 7` (Git Bash on Windows, CuPy and an NVIDIA GPU; `kformula.py` imports the A000944 formula modules from `../A000944/formula`, published alongside): writes the cycle-index rows
   `kf_30_<prime>.json` for 7 primes below 2^24 (about 13 min per prime on an RTX 3080), then runs
   `tower_gpu.py 30 7`, which prints the three sequences with CRT and writes `tower_30.json`.
   Without a GPU: `python kformula.py 30 <prime_index>` for each prime, then `python run_tower.py 30 7`
   (slow: hours per prime).
2. Checks: `python labeled_check.py 30 7`, `python marginals_check.py 30 7`, and in `brute/` the plantri
   plug-in (`zkplug.c`, `zkplug_compare.py`).

## What is being counted
A005470(n) is the number of simple planar graphs on n unlabeled vertices; A003094 counts the connected
ones and A021103 the 2-connected ones (blocks). For n = 4 the 11 graphs are all graphs on 4 vertices.

## Method: the planar-graph tower
Every step but the first is classical species theory. The first step uses the A000944 formula machinery,
so nothing is generated.

1. **3-connected planar graphs, full cycle index** (`kformula.py`). By Whitney's theorem their automorphisms
   are the map automorphisms of the unique embedding. So the two-sort cycle index
   Z_K = Σ_K 1/|Aut K| Σ_g (vertex cycles of g)·(edge orbits of g, each marked "flipped" when g^k reverses
   its edges) is a sum over the geometric types of A000944 (`../A000944/README.md`). The cycle structure of
   g is forced by its type, except that reflections (and C_mh) need the numbers of mirror vertices, mirror
   edges and crossing edges separately. The half-map grammar (`disk.py`) is therefore run with three
   markers: on V labels, on X labels, and on V→V boundary edges (= mirror edges). It is evaluated on a
   grid and interpolated. Per type:

   | type | vertex cycles | edge orbits |
   |---|---|---|
   | identity | p1^v | (1,0)^e |
   | rotation k ≥ 3 | p1^#Vpoles · pk^rest | (k,0)^(e/k) |
   | half-turn | p1^#V p2^rest | (1,1)^#E (2,0)^rest |
   | reflection | p1^nv p2^rest | (1,0)^nm (1,1)^nc (2,0)^rest |
   | antipodal | p2^(v/2) | (2,0)^(e/2) |
   | rot. refl. order n | p2^[VV] pn^rest | (2,1)^[EE] (n,0)^rest |
   | C_mh (n = 2m) | pm^(nv/m) p2^[V poles] p2m^rest | (m,0)^(nm/m) (m,1)^(nc/m) (2m,0)^rest |

2. **Networks** (`tower_p.py`). Two-pole networks are Z2-species (pole swap). Each class carries a
   pole-fixing cycle index F^id and a pole-swapping one F^τ: D = e + S + P + H, where S is a series chain of
   ≥ 2 networks from D − S (reversal for τ), P is a set of ≥ 2 from S + H plus at most one edge, and H is a
   3-connected core with its root edge deleted and every other edge replaced by D. A flipped edge orbit of
   length k takes D^τ[k] (plethysm p_j → p_jk), an unflipped one D^id[k]. H^id and H^τ are the edge
   derivatives of Z_K. The substitution is expanded in binomial moments (D = 1 + d with d of weight ≥ 1),
   so only low powers of d enter.
3. **2-connected graphs** by the dissymmetry theorem on the SPQR tree:
   Z_B = K2 + R + S-node + P-node − SP − SR − PR + (RR unordered − RR ordered), with dihedral cycle
   indices for the S-nodes and bond (pole-swap) indices for the P-nodes.
4. **Connected graphs** (Robinson / Harary–Palmer): C• = x·exp(Σ_k B'(C•)(x^k)/k) and
   C = C• + B(C•) − C•·B'(C•). **All graphs**: Euler transform.

All arithmetic is mod primes below 2^24 with CRT; the tower works with dense vectors over the
partitions of weight ≤ N.

## Checks
* `brute/zk.py`: the cycle index of all 3-connected planar graphs with n ≤ 11, by brute force (nauty
  graphs from `../A000944/brute`, every automorphism enumerated). Its vertex sums are A000944(n).
* `tower.py` (exact Fractions) and `tower_p.py` (mod p) on the brute-force cycle index reproduce
  A021103, A003094 and A005470 exactly for n ≤ 11.
* `kformula.py check 11`: the formula cycle index equals the brute-force one on all 637 cells
  (n, vertex cycle type, edge orbit type) with n ≤ 11, all symmetry types included.
* `labeled_check.py 30 7` (identity sector, mod each of the 7 primes): the p1^n coefficients of the
  formula cycle index and of the tower's block series B give n!-scaled labeled counts, and a labeled
  Robinson/exp step gives the rest. All agree with the published b-files (Bodirsky, Gröpl, Kang) for
  every n ≤ 30: A096330 (3-connected), A096331 (2-connected), A096332 (connected), A066537 (all).
  This verifies the networks, SPQR dissymmetry and identity rows at full size; the projection p_k -> 0
  (k ≥ 2) kills every non-identity term, so the symmetric sectors are not tested by it.
* `marginals_check.py 30 7`: summing the formula rows over cycle types reproduces exhaustive-generation
  (plantri) counts mod each prime: A000944 for n ≤ 18, A002840 by edges for e ≤ 35, the whole
  A049337 (n, e) triangle for n ≤ 14 (462 cells), A000109 triangulations for n ≤ 24. This checks
  every symmetry type's total weight beyond n = 11, though not its split by cycle type.
* `brute/zkplug.c`: a plantri 5.8 plug-in (`cc -o plantri_zk -O3 -DPLUGIN=zkplug.c plantri.c`, run as
  `plantri_zk -pG n res/mod`; `brute/plug/run_split.sh n parts out` runs 6 parts in WSL) that reads the
  full automorphism group of every polyhedron from plantri and records each automorphism's vertex
  cycle type and edge orbits. It equals `brute/zk.py` exactly for n ≤ 11, and
  `brute/zkplug_compare.py n plug/zk_n.txt --kf 30 7` finds the formula rows equal to it on every cell,
  mod all 7 primes, for n = 12, 13, 14 and 15 (322, 371, 556 and 673 cells; n = 15 is 23.8 billion
  polyhedra, 30 min on 6 cores). So the 3-connected cycle index is verified cell by cell through n = 15.

## Files
* `brute/zk.py`, `brute/zk_le10.json`, `brute/zk_11.json`: brute-force cycle index, n ≤ 11.
* `tower.py`: exact pilot of the tower; `tower_p.py`: the fast mod-p tower.
* `kformula.py`: the formula cycle index (`python kformula.py N prime_index` writes `kf_N_p.json`).
* `run_tower.py N nprimes`: tower on the formula rows, CRT, comparison with the known terms.
* `tower_gpu.py N nprimes`: the same tower on the GPU (segmented-sum product kernel, moment tables as
  row-sum + scatter kernels, vectorized univariate stage); same output as `run_tower.py`.
* `gpu.py`: CuPy series class for the reflection grammar. Its multiply is a direct, exact truncated 2D
  convolution in a CUDA kernel, one block per grid point (`python gpu.py` checks it against fastmul).
  Enable with `KF_GPU=1` (optional `KF_GPU_CHUNK`, default 1024 points per batch; `KF_OUT` sets the
  output folder). `run_gpu.sh N NPRIMES`: the GPU run (two streams of primes, then the tower).

## Log
- 2026-11-05: pilot tower validated on brute-force 3-connected data (n ≤ 11); formula cycle index
  validated cell by cell (n ≤ 11); N = 15 formula run started.
- 2026-11-05: speedups. `ggram.py` is a dedicated GPU port of the half-map grammar: Barrett reduction
  instead of int64 `%` (which costs more than the convolution on a GPU), monomials as shift kernels,
  host-side lower degree bounds (no syncs), and pass t of the arc fixed point computed only up to
  x-degree t. The marker grid is a lower set (nv + 2nc ≤ VX, nm + 2nc ≤ VX) with Newton interpolation,
  1771 points instead of 4851 at N = 20. The tower solves the networks weight by weight. All of this
  gives rows/terms identical to the earlier runs at N = 15 and 20. Per prime: N = 20 takes about 1 min
  for the rows (the first GPU version took 25 min, the CPU about 2.2 h); N = 25 takes 3–5 min for the rows
  and about 6 min for the tower (`run_gpu.sh N NPRIMES`, `KF_GPU=2`).
- 2026-11-05: **N = 25** (`tower_25.log`, 5 primes, CRT stable with 4). It agrees with the N = 20 run on
  every common term. A005470(21..25) = 266892197689888343244, 6084033586294190154202, 140012314430574390754809, 3249784959167329158913333, 76016793864553053210040815;
  A003094(21..25) = 254499290386523690634, 5804661393043936816211, 133647000805745363204565, 3103364947739587802644026, 72619635686723255869114898;
  A021103(21..25) = 96009973957009506442, 2109220961695218583498, 46763841586042636362681, 1045438229538044361911334, 23548181230365705962266586.
  a(n)/a(n-1)·(n/(n-1))^(7/2) for A005470 rises smoothly: 26.11 (n = 16) to 26.98 (n = 25).
- 2026-11-05: tower on the GPU (`tower_gpu.py`): N = 20 takes 12.5 s for 5 primes (CPU 69 s), N = 25 takes
  30 s (CPU about 30 min while the machine was shared). `tower_20.json` and `tower_25.json` are
  byte-identical to the CPU results (`tg25.log`).
- 2026-11-06: **N = 30** (`tower_30.log`, 7 primes, CRT stable with 6; values up to 2^109). 68206 cycle-index
  rows per prime, about 13 min per prime with two streams on the GPU (`KF_GPU_CHUNK=512`), 50 min in all;
  the GPU tower takes seconds. It agrees with the N = 25 run on every common term.
  A005470(26..30) = 1790733466003422184670125640, 42457813673889602379288800530, 1012653135463833979592320804807, 24285021128825709709135132968024, 585345393678411726175403670215562;
  A003094(26..30) = 1711297840974659008783320156, 40587146973195981103012285191, 968313129216044117167721499740, 23227758888209768266519857560467, 559996781949783047710910447811156;
  A021103(26..30) = 534073930971881285532076914, 12189355449099428544303427991, 279818237785356119481685635143, 6457929182734030291113816963729, 149782224970101334436971067029495.
  The corrected growth ratio of A005470 reaches 27.14 at n = 30, still rising smoothly.
- 2026-11-06: independent checks (`labeled_check.py`, `marginals_check.py`, see Checks): all labeled counts
  to n = 30 and all plantri marginals agree. Then the plantri plug-in brute force (`brute/zkplug.c`):
  the formula cycle index matches cell by cell for every n ≤ 15. What remains unverified independently
  is only the cycle-type split of the symmetric sectors for 16 ≤ n ≤ 30, which the formula computes
  with the same code that is exact through n = 15.
