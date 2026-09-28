//! Independent sets of the hypercube Q_n (OEIS A027624), via Q_{d+2} = Q_d x C_4.
//!
//! An independent set of Q_{d+2} is a ring S_0, S_1, S_2, S_3 of independent sets of Q_d with
//! neighbouring members disjoint. Summing out S_1 and S_3:
//!
//!     a(d+2) = sum over independent S, U of Q_d of  f(~(S | U))^2,
//!
//! where f(W) = number of independent sets of Q_d inside the vertex set W. Likewise
//! Q_{d+1} = Q_d x K_2 gives a(d+1) = sum_S f(~S), which is used as a check.
//!
//! f is evaluated through Q_d = Q_{d-1} x K_2: an independent set of Q_d is a pair (A, B) of
//! independent sets of Q_{d-1} with A, B disjoint, so
//!
//!     f(W) = sum over independent A of Q_{d-1} with A inside W_lo of  g(W_hi & ~A),
//!
//! with g the subset-sum ("zeta") table of Q_{d-1} (2^16 entries for d = 5) and, for every W_lo,
//! a precomputed list of the independent sets inside it.
//!
//! The summand is invariant under Aut(Q_d) acting on S and U together, so S runs over orbit
//! representatives only, weighted by orbit size.
//!
//! Usage:
//!   hypercube-indep validate                       d = 1..4, every group, against known terms
//!   hypercube-indep run <d> [--group full|flips|none] [--threads T] [--progress SECS]
//!   hypercube-indep spotcheck <d> [samples]        f against a direct count on random W
//!   hypercube-indep a8 ... | a8run [--nosym]        a(8) on the GPU, natively (see gpu.rs)
//!   hypercube-indep poly <d> [--group ...]         same sums by set size: the independence
//!                                                  polynomials of Q_{d+1} and Q_{d+2}
//!                                                  (rows of A354802; value at -1 is A354082)

mod gpu;
mod ygroup;

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::{Duration, Instant};

/// a(0..6) from the OEIS entry.
const KNOWN: [u128; 7] = [2, 3, 7, 35, 743, 254475, 19768832143];

/// A354082(0..7): the independence polynomial of Q_n at x = -1.
const KNOWN_ALT: [i128; 8] = [0, -1, -1, 3, 7, 11, 143, 7715];

/// A354802 DATA (rows n = 0..5 in full, then the start of row 6): independent sets by size.
const KNOWN_ROWS: &[u64] = &[
    1, 1, 1, 2, 1, 4, 2, 1, 8, 16, 8, 2, 1, 16, 88, 208, 228, 128, 56, 16, 2, 1, 32, 416, 2880, 11760,
    29856, 48960, 54304, 44240, 29920, 17952, 9088, 3672, 1120, 240, 32, 2, 1, 64, 1824, 30720, 342400,
    2682432, 15328064, 65515840, 213464240, 538811200, 1070860384, 1708551424, 2245780976, 2517976640,
    2509047680,
];

// ---------------------------------------------------------------------------------------------
// Independent sets and the f evaluator
// ---------------------------------------------------------------------------------------------

fn full_mask(nv: usize) -> u32 {
    if nv == 32 { u32::MAX } else { (1u32 << nv) - 1 }
}

/// All independent sets of Q_d as vertex bitmasks (vertex v = bit v), sorted.
fn indep_sets(d: usize) -> Vec<u32> {
    let nv = 1usize << d;
    let nbr: Vec<u32> = (0..nv).map(|v| (0..d).fold(0u32, |m, i| m | 1 << (v ^ (1 << i)))).collect();
    fn rec(v: usize, cur: u32, nbr: &[u32], out: &mut Vec<u32>) {
        if v == nbr.len() {
            out.push(cur);
            return;
        }
        rec(v + 1, cur, nbr, out);
        if cur & nbr[v] == 0 {
            rec(v + 1, cur | 1 << v, nbr, out);
        }
    }
    let mut out = vec![];
    rec(0, 0, &nbr, &mut out);
    out.sort_unstable();
    out
}

/// Counts independent sets of Q_d inside W through the split Q_d = Q_{d-1} x K_2.
struct FEval {
    half: usize,     // vertices per half = 2^(d-1)
    hmask: u32,      // mask of one half
    g: Vec<u32>,     // g[X] = # independent sets of Q_{d-1} inside X
    off: Vec<u32>,   // CSR offsets over X_lo
    lst: Vec<u32>,   // independent sets of Q_{d-1} inside X_lo
}

impl FEval {
    fn new(d: usize) -> FEval {
        assert!((1..=5).contains(&d));
        let half = 1usize << (d - 1);
        let hmask = full_mask(half);
        let sub = indep_sets(d - 1);
        let size = 1usize << half;
        let mut g = vec![0u32; size];
        for &a in &sub {
            g[a as usize] = 1;
        }
        for b in 0..half {
            for x in 0..size {
                if x >> b & 1 == 1 {
                    g[x] += g[x ^ 1 << b];
                }
            }
        }
        // CSR: for every X, the independent sets A with A inside X
        let mut cnt = vec![0u32; size + 1];
        for &a in &sub {
            let comp = !a & hmask;
            let mut s = comp;
            loop {
                cnt[(a | s) as usize + 1] += 1;
                if s == 0 {
                    break;
                }
                s = (s - 1) & comp;
            }
        }
        for x in 0..size {
            cnt[x + 1] += cnt[x];
        }
        let mut fill = cnt.clone();
        let mut lst = vec![0u32; cnt[size] as usize];
        for &a in &sub {
            let comp = !a & hmask;
            let mut s = comp;
            loop {
                let x = (a | s) as usize;
                lst[fill[x] as usize] = a;
                fill[x] += 1;
                if s == 0 {
                    break;
                }
                s = (s - 1) & comp;
            }
        }
        FEval { half, hmask, g, off: cnt, lst }
    }

    #[inline]
    fn f(&self, w: u32) -> u64 {
        let lo = (w & self.hmask) as usize;
        let hi = (w >> self.half) & self.hmask;
        let (a, b) = (self.off[lo] as usize, self.off[lo + 1] as usize);
        let mut s = 0u64;
        for &x in &self.lst[a..b] {
            s += self.g[(hi & !x) as usize] as u64;
        }
        s
    }
}

/// Direct count of independent sets of Q_d inside W, for spot checks.
fn f_direct(d: usize, w: u32) -> u64 {
    let nv = 1usize << d;
    let verts: Vec<usize> = (0..nv).filter(|&v| w >> v & 1 == 1).collect();
    fn rec(i: usize, cur: u32, verts: &[usize], d: usize) -> u64 {
        if i == verts.len() {
            return 1;
        }
        let v = verts[i];
        let mut n = rec(i + 1, cur, verts, d);
        if (0..d).all(|b| cur >> (v ^ 1 << b) & 1 == 0) {
            n += rec(i + 1, cur | 1 << v, verts, d);
        }
        n
    }
    rec(0, 0, &verts, d)
}

// ---------------------------------------------------------------------------------------------
// Symmetry: orbits of Aut(Q_d) (or a subgroup) on independent sets
// ---------------------------------------------------------------------------------------------

/// A vertex permutation of Q_d, applied to bitmasks through per-byte lookup tables.
pub(crate) struct VPerm {
    tab: [[u32; 256]; 4],
}

impl VPerm {
    fn new(p: &[usize]) -> VPerm {
        let mut tab = [[0u32; 256]; 4];
        for (k, t) in tab.iter_mut().enumerate() {
            for (byte, e) in t.iter_mut().enumerate() {
                for bit in 0..8 {
                    let v = 8 * k + bit;
                    if byte >> bit & 1 == 1 && v < p.len() {
                        *e |= 1 << p[v];
                    }
                }
            }
        }
        VPerm { tab }
    }
    #[inline]
    fn apply(&self, m: u32) -> u32 {
        self.tab[0][(m & 255) as usize]
            | self.tab[1][(m >> 8 & 255) as usize]
            | self.tab[2][(m >> 16 & 255) as usize]
            | self.tab[3][(m >> 24) as usize]
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Group {
    Full,  // Aut(Q_d): bit flips and coordinate permutations, order 2^d d!
    Flips, // bit flips only, order 2^d
    None,  // trivial group
}

impl Group {
    fn parse(s: &str) -> Group {
        match s {
            "full" => Group::Full,
            "flips" => Group::Flips,
            "none" => Group::None,
            _ => panic!("unknown group {s}"),
        }
    }
}

/// Generators of the chosen group as vertex permutations of Q_d.
fn generators(d: usize, grp: Group) -> Vec<VPerm> {
    let nv = 1usize << d;
    let mut gens = vec![];
    if grp == Group::None {
        return gens;
    }
    for i in 0..d {
        let p: Vec<usize> = (0..nv).map(|v| v ^ 1 << i).collect();
        gens.push(VPerm::new(&p));
    }
    if grp == Group::Full && d >= 2 {
        // coordinate permutations are generated by the transposition (0 1) and the d-cycle
        let permute = |sigma: &dyn Fn(usize) -> usize| -> Vec<usize> {
            (0..nv).map(|v| (0..d).fold(0, |w, i| w | (v >> i & 1) << sigma(i))).collect()
        };
        gens.push(VPerm::new(&permute(&|i| if i < 2 { 1 - i } else { i })));
        gens.push(VPerm::new(&permute(&|i| (i + 1) % d)));
    }
    gens
}

/// Orbit representatives (indices into `sets`) and orbit sizes, by union-find over generators.
fn orbits(sets: &[u32], gens: &[VPerm]) -> Vec<(usize, u64)> {
    let n = sets.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for (i, &s) in sets.iter().enumerate() {
        for g in gens {
            let j = sets.binary_search(&g.apply(s)).expect("image of an independent set is independent");
            let (a, b) = (find(&mut parent, i), find(&mut parent, j));
            if a != b {
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut size = vec![0u64; n];
    for i in 0..n {
        let r = find(&mut parent, i);
        size[r] += 1;
    }
    (0..n).filter(|&i| parent[i] == i).map(|i| (i, size[i])).collect()
}

// ---------------------------------------------------------------------------------------------
// The two sums
// ---------------------------------------------------------------------------------------------

static P_DONE: AtomicUsize = AtomicUsize::new(0);
static P_TOTAL: AtomicUsize = AtomicUsize::new(0);
static P_STOP: AtomicBool = AtomicBool::new(false);
static THREADS: AtomicUsize = AtomicUsize::new(0);

pub(crate) fn num_threads() -> usize {
    match THREADS.load(Ordering::Relaxed) {
        // leave two logical CPUs for the user's other jobs
        0 => std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).saturating_sub(2).max(1),
        t => t,
    }
}

fn spawn_heartbeat(secs: u64) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        let t0 = Instant::now();
        let mut next = Duration::from_secs(secs);
        while !P_STOP.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_millis(200));
            if t0.elapsed() < next {
                continue;
            }
            next += Duration::from_secs(secs);
            let (done, tot) = (P_DONE.load(Ordering::Relaxed), P_TOTAL.load(Ordering::Relaxed));
            if tot == 0 {
                continue;
            }
            let el = t0.elapsed().as_secs_f64();
            let eta = if done > 0 { format!("~{:.0}s", el / done as f64 * (tot - done) as f64) } else { "?".into() };
            eprintln!("[{:>6.0}s] items {done}/{tot} ({:.1}%)  ETA {eta}", el, 100.0 * done as f64 / tot as f64);
        }
    })
}

/// a(d+1) = sum_S f(~S)   (Q_{d+1} = Q_d x K_2).
fn sum_k2(fe: &FEval, sets: &[u32], full: u32) -> u128 {
    sets.iter().map(|&s| fe.f(!s & full) as u128).sum()
}

/// a(d+2) = sum over orbit reps S (weighted) and all U of f(~(S|U))^2   (Q_{d+2} = Q_d x C_4).
fn sum_c4(fe: &FEval, sets: &[u32], reps: &[(usize, u64)], full: u32) -> u128 {
    const CHUNK: usize = 4096;
    let nchunks = sets.len().div_ceil(CHUNK);
    let nitems = reps.len() * nchunks;
    P_DONE.store(0, Ordering::Relaxed);
    P_TOTAL.store(nitems, Ordering::Relaxed);
    let next = AtomicUsize::new(0);
    let total = Mutex::new(0u128);
    std::thread::scope(|scope| {
        for _ in 0..num_threads() {
            scope.spawn(|| {
                let mut local = 0u128;
                loop {
                    let it = next.fetch_add(1, Ordering::Relaxed);
                    if it >= nitems {
                        break;
                    }
                    let (ri, ci) = (it / nchunks, it % nchunks);
                    let (si, w) = reps[ri];
                    let s = sets[si];
                    let mut inner = 0u64; // <= CHUNK * 254475^2 < 2^49
                    for &u in &sets[ci * CHUNK..((ci + 1) * CHUNK).min(sets.len())] {
                        let x = fe.f(!(s | u) & full);
                        inner += x * x;
                    }
                    local += w as u128 * inner as u128;
                    P_DONE.fetch_add(1, Ordering::Relaxed);
                }
                *total.lock().unwrap() += local;
            });
        }
    });
    total.into_inner().unwrap()
}

struct RunResult {
    k2: u128,
    c4: u128,
    nsets: usize,
    nreps: usize,
    csr_len: usize,
    secs: f64,
}

fn run(d: usize, grp: Group, verbose: bool) -> RunResult {
    let t0 = Instant::now();
    let nv = 1usize << d;
    let full = full_mask(nv);
    let sets = indep_sets(d);
    let fe = FEval::new(d);
    let gens = generators(d, grp);
    let reps = orbits(&sets, &gens);
    let wsum: u64 = reps.iter().map(|r| r.1).sum();
    assert_eq!(wsum as usize, sets.len(), "orbit sizes must add up to the number of sets");
    if verbose {
        eprintln!(
            "d={d} group={grp:?}: {} independent sets, {} orbit reps, CSR {} entries, setup {:.2}s",
            sets.len(),
            reps.len(),
            fe.lst.len(),
            t0.elapsed().as_secs_f64()
        );
    }
    let k2 = sum_k2(&fe, &sets, full);
    let c4 = sum_c4(&fe, &sets, &reps, full);
    RunResult { k2, c4, nsets: sets.len(), nreps: reps.len(), csr_len: fe.lst.len(), secs: t0.elapsed().as_secs_f64() }
}

fn check(label: &str, got: u128, n: usize) -> bool {
    match KNOWN.get(n) {
        Some(&k) => {
            let ok = k == got;
            println!("  {label}: a({n}) = {got}  {}", if ok { "OK" } else { "MISMATCH" });
            ok
        }
        None => {
            println!("  {label}: a({n}) = {got}  (NEW)");
            true
        }
    }
}

// ---------------------------------------------------------------------------------------------
// The same sums refined by set size (independence polynomials)
// ---------------------------------------------------------------------------------------------

const GC: usize = 9; // coefficients of G: independent sets of Q_{d-1} have size <= 8 (d <= 5)
const FC: usize = 17; // coefficients of F: size <= 16 in Q_d
const PC: usize = 65; // coefficients of the Q_{d+2} polynomial: size <= 64

/// F(W) = sum over independent T of Q_d inside W of x^|T|, through Q_d = Q_{d-1} x K_2.
struct FPoly {
    half: usize,
    hmask: u32,
    g: Vec<[u32; GC]>, // G[X] = sum over independent A of Q_{d-1} inside X of x^|A|
    off: Vec<u32>,     // same CSR lists as FEval
    lst: Vec<u32>,
}

impl FPoly {
    fn new(d: usize) -> FPoly {
        let fe = FEval::new(d);
        let size = 1usize << fe.half;
        let mut g = vec![[0u32; GC]; size];
        for &a in &indep_sets(d - 1) {
            g[a as usize][a.count_ones() as usize] = 1;
        }
        for b in 0..fe.half {
            for x in 0..size {
                if x >> b & 1 == 1 {
                    let lower = g[x ^ 1 << b];
                    for k in 0..GC {
                        g[x][k] += lower[k];
                    }
                }
            }
        }
        FPoly { half: fe.half, hmask: fe.hmask, g, off: fe.off, lst: fe.lst }
    }

    #[inline]
    fn f(&self, w: u32) -> [u64; FC] {
        let lo = (w & self.hmask) as usize;
        let hi = (w >> self.half) & self.hmask;
        let mut out = [0u64; FC];
        for &a in &self.lst[self.off[lo] as usize..self.off[lo + 1] as usize] {
            let k = a.count_ones() as usize;
            let row = &self.g[(hi & !a) as usize];
            for j in 0..GC {
                out[k + j] += row[j] as u64;
            }
        }
        out
    }
}

/// P_{d+1}(x) = sum_S x^|S| F(~S)   (Q_{d+1} = Q_d x K_2).
fn poly_k2(fp: &FPoly, sets: &[u32], full: u32) -> Vec<u128> {
    let mut p = vec![0u128; PC];
    for &s in sets {
        let k = s.count_ones() as usize;
        for (j, c) in fp.f(!s & full).iter().enumerate() {
            p[k + j] += *c as u128;
        }
    }
    p
}

/// P_{d+2}(x) = sum over orbit reps S (weighted) and all U of x^(|S|+|U|) F(~(S|U))^2.
fn poly_c4(fp: &FPoly, sets: &[u32], reps: &[(usize, u64)], full: u32) -> Vec<u128> {
    const CHUNK: usize = 4096;
    let nchunks = sets.len().div_ceil(CHUNK);
    let nitems = reps.len() * nchunks;
    P_DONE.store(0, Ordering::Relaxed);
    P_TOTAL.store(nitems, Ordering::Relaxed);
    let next = AtomicUsize::new(0);
    let total = Mutex::new(vec![0u128; PC]);
    std::thread::scope(|scope| {
        for _ in 0..num_threads() {
            scope.spawn(|| {
                let mut local = vec![0u128; PC];
                loop {
                    let it = next.fetch_add(1, Ordering::Relaxed);
                    if it >= nitems {
                        break;
                    }
                    let (ri, ci) = (it / nchunks, it % nchunks);
                    let (si, w) = reps[ri];
                    let s = sets[si];
                    let mut inner = [0u64; PC]; // coefficients <= CHUNK * f^2 < 2^49
                    for &u in &sets[ci * CHUNK..((ci + 1) * CHUNK).min(sets.len())] {
                        let f = fp.f(!(s | u) & full);
                        let sh = (s.count_ones() + u.count_ones()) as usize;
                        for i in 0..FC {
                            if f[i] == 0 {
                                continue;
                            }
                            for j in 0..FC {
                                inner[sh + i + j] += f[i] * f[j];
                            }
                        }
                    }
                    for k in 0..PC {
                        local[k] += w as u128 * inner[k] as u128;
                    }
                    P_DONE.fetch_add(1, Ordering::Relaxed);
                }
                let mut t = total.lock().unwrap();
                for k in 0..PC {
                    t[k] += local[k];
                }
            });
        }
    });
    total.into_inner().unwrap()
}

/// Drops trailing zero coefficients.
fn trim(p: &[u128]) -> Vec<u128> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}

/// Row n of A354802 from KNOWN_ROWS, if it is there in full (row n has 2^(n-1) + 1 entries).
fn known_row(n: usize) -> Option<Vec<u128>> {
    let len = |m: usize| if m == 0 { 2 } else { (1usize << (m - 1)) + 1 };
    let start: usize = (0..n).map(len).sum();
    let end = start + len(n);
    if end <= KNOWN_ROWS.len() {
        Some(KNOWN_ROWS[start..end].iter().map(|&x| x as u128).collect())
    } else {
        None
    }
}

/// Checks the polynomial of Q_n against A027624 (sum), A354082 (value at -1) and A354802 (row).
fn check_poly(label: &str, p: &[u128], n: usize) -> bool {
    let p = trim(p);
    let total: u128 = p.iter().sum();
    let alt: i128 = p.iter().enumerate().map(|(k, &c)| if k % 2 == 0 { c as i128 } else { -(c as i128) }).sum();
    let mut ok = true;
    let mut notes = vec![];
    match KNOWN.get(n) {
        Some(&k) => {
            ok &= k == total;
            notes.push(if k == total { "sum OK" } else { "sum MISMATCH" });
        }
        None => notes.push("sum NEW"),
    }
    if let Some(&k) = KNOWN_ALT.get(n) {
        ok &= k == alt;
        notes.push(if k == alt { "A354082 OK" } else { "A354082 MISMATCH" });
    }
    match known_row(n) {
        Some(r) => {
            ok &= r == p;
            notes.push(if r == p { "A354802 row OK" } else { "A354802 row MISMATCH" });
        }
        None => notes.push("A354802 row not in DATA"),
    }
    println!("  {label}: n={n}  sum={total}  P(-1)={alt}  [{}]", notes.join(", "));
    println!("    row: {}", p.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","));
    ok
}

fn run_poly(d: usize, grp: Group) -> bool {
    let t0 = Instant::now();
    let full = full_mask(1 << d);
    let sets = indep_sets(d);
    let fp = FPoly::new(d);
    let reps = orbits(&sets, &generators(d, grp));
    let k2 = poly_k2(&fp, &sets, full);
    let c4 = poly_c4(&fp, &sets, &reps, full);
    println!("poly d={d} group={grp:?}: {} reps, {:.2}s", reps.len(), t0.elapsed().as_secs_f64());
    let a = check_poly("Q_d x K_2", &k2, d + 1);
    let b = check_poly("Q_d x C_4", &c4, d + 2);
    a && b
}

// ---------------------------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------------------------

fn validate() -> bool {
    let mut ok = true;
    for d in 1..=4 {
        for grp in [Group::Full, Group::Flips, Group::None] {
            let r = run(d, grp, false);
            println!("d={d} group={grp:?} ({} sets, {} reps, {:.3}s)", r.nsets, r.nreps, r.secs);
            ok &= check("Q_d x K_2", r.k2, d + 1);
            ok &= check("Q_d x C_4", r.c4, d + 2);
            ok &= KNOWN[d] == r.nsets as u128;
        }
    }
    println!("{}", if ok { "ALL OK" } else { "FAILURES" });
    ok
}

fn spotcheck(d: usize, samples: usize) -> bool {
    let fe = FEval::new(d);
    let sets = indep_sets(d);
    let full = full_mask(1 << d);
    let mut rng = 0x9E3779B97F4A7C15u64;
    let mut next = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    let mut bad = 0;
    for i in 0..samples {
        // half the samples are the W that the C_4 sum actually meets, half are uniform
        let w = if i % 2 == 0 {
            let s = sets[next() as usize % sets.len()];
            let u = sets[next() as usize % sets.len()];
            !(s | u) & full
        } else {
            next() as u32 & full
        };
        let (a, b) = (fe.f(w), f_direct(d, w));
        if a != b {
            bad += 1;
            eprintln!("mismatch W={w:#010x}: csr {a}, direct {b}");
        }
    }
    println!("spotcheck d={d}: {samples} samples, {bad} mismatches");
    bad == 0
}

/// Windows 11 treats console jobs started in the background as "efficiency" work (EcoQoS) and keeps
/// them on the E-cores of hybrid CPUs (e.g. i5-12600K: 6 P-cores idle, 4 E-cores at 100%).
/// Opt this process out of execution-speed throttling so the scheduler may use the P-cores.
#[cfg(windows)]
fn disable_power_throttling() {
    #[repr(C)]
    struct ProcessPowerThrottlingState {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
    }
    const PROCESS_POWER_THROTTLING: i32 = 4; // PROCESS_INFORMATION_CLASS::ProcessPowerThrottling
    const EXECUTION_SPEED: u32 = 0x1; // PROCESS_POWER_THROTTLING_EXECUTION_SPEED
    let st = ProcessPowerThrottlingState { version: 1, control_mask: EXECUTION_SPEED, state_mask: 0 };
    let ok = unsafe {
        SetProcessInformation(
            GetCurrentProcess(),
            PROCESS_POWER_THROTTLING,
            &st as *const _ as *const std::ffi::c_void,
            std::mem::size_of::<ProcessPowerThrottlingState>() as u32,
        )
    };
    if ok == 0 {
        eprintln!("note: could not disable Windows power throttling (EcoQoS)");
    }
}

#[cfg(not(windows))]
fn disable_power_throttling() {}

fn main() {
    disable_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let opt = |name: &str| args.iter().position(|a| a == name).map(|i| args[i + 1].clone());
    if let Some(t) = opt("--threads") {
        THREADS.store(t.parse().unwrap(), Ordering::Relaxed);
    }
    let every: u64 = opt("--progress").map(|s| s.parse().unwrap()).unwrap_or(30);
    let ok = match args.first().map(String::as_str) {
        Some("validate") => validate(),
        Some("spotcheck") => {
            let d: usize = args[1].parse().unwrap();
            let n: usize = args.get(2).filter(|a| !a.starts_with("--")).map(|a| a.parse().unwrap()).unwrap_or(2000);
            spotcheck(d, n)
        }
        Some("a8") => gpu::a8(&args[1..]),
        Some("a8run") => gpu::a8run(&args[1..]),
        Some("reps6") => {
            let path = std::path::PathBuf::from(opt("--out").unwrap_or_else(|| "results/ylo_reps_d6.bin".into()));
            std::fs::create_dir_all(path.parent().unwrap()).ok();
            let reps = ygroup::enumerate_reps6(&path);
            reps.len() == 1_228_158
        }
        Some("cpu6") => {
            // cpu6 <rep index> ... : reference rep sums for GPU validation
            let reps = ygroup::load_reps6(std::path::Path::new(&opt("--reps").unwrap_or_else(|| "results/ylo_reps_d6.bin".into())));
            let idx: Vec<usize> = args[1..].iter().take_while(|a| !a.starts_with("--")).map(|a| a.parse().unwrap()).collect();
            let sums = ygroup::cpu_reps6(&reps, &idx, args.iter().any(|a| a == "--sym"));
            for (i, s) in idx.iter().zip(&sums) {
                println!("{} {} {} {}", i, reps[*i].0, reps[*i].1, s.to_decimal());
            }
            true
        }
        Some("ybench") => {
            let n: usize = args.get(1).map(|a| a.parse().unwrap()).unwrap_or(8);
            let y3: u64 = args.get(2).map(|a| a.parse().unwrap()).unwrap_or(8);
            ygroup::bench6(n, y3);
            true
        }
        Some("ygroup") => {
            // union-grouped sum, d = 3..5 -> a(5..7)
            let mut ok = true;
            for d in 3..=5 {
                let t = ygroup::run(d);
                let exp = KNOWN.get(d + 2).copied();
                let got = t.to_decimal();
                let good = exp.map(|e| ygroup::U256 { hi: 0, lo: e } == t).unwrap_or(true);
                println!("ygroup d={d}: a({}) = {got}  {}", d + 2, if good { "OK" } else { "MISMATCH" });
                ok &= good;
            }
            ok
        }
        Some("poly") => {
            let d: usize = args[1].parse().unwrap();
            let grp = Group::parse(&opt("--group").unwrap_or_else(|| "full".into()));
            let hb = if every > 0 { Some(spawn_heartbeat(every)) } else { None };
            let ok = if d == 0 { (1..=5).all(|d| run_poly(d, grp)) } else { run_poly(d, grp) };
            P_STOP.store(true, Ordering::Relaxed);
            if let Some(h) = hb {
                h.join().ok();
            }
            ok
        }
        Some("run") => {
            let d: usize = args[1].parse().unwrap();
            let grp = Group::parse(&opt("--group").unwrap_or_else(|| "full".into()));
            eprintln!("threads: {}", num_threads());
            let hb = if every > 0 { Some(spawn_heartbeat(every)) } else { None };
            let r = run(d, grp, true);
            P_STOP.store(true, Ordering::Relaxed);
            if let Some(h) = hb {
                h.join().ok();
            }
            println!("d={d} group={grp:?}: {} sets, {} reps, CSR {} entries, {:.2}s", r.nsets, r.nreps, r.csr_len, r.secs);
            let a = check("Q_d x K_2", r.k2, d + 1);
            let b = check("Q_d x C_4", r.c4, d + 2);
            a && b
        }
        _ => {
            eprintln!("usage: hypercube-indep validate | run <d> [--group full|flips|none] | poly <d, or 0 for 1..5> | spotcheck <d> [samples]");
            false
        }
    };
    std::process::exit(if ok { 0 } else { 1 });
}
