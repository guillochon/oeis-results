//! Burnside engine for OEIS A109455: threshold functions of n variables up to permuting variables.
//! Compiled port of ../engine.py (see that file and ../README.md for the mathematics).
//!
//!   a(n) = (1/n!) * sum over cycle types lam of |class(lam)| * Fix(lam)
//!
//! Fix(lam) = number of threshold functions on the grid {0..l_1} x ... x {0..l_k} (l_j = cycle
//! lengths), counted as sum over *regular positive* threshold functions f of
//! orbit(f) * 2^(#axes f depends on). Regular positive functions are enumerated point by point in a
//! linear extension of the dominance order; a carried witness (w, t) decides one branch for free and
//! an LP is solved only for the other one.
//!
//! LP: we need (w, t) with  w.x >= t on 1-points, w.x <= t - 1 on 0-points, w >= 0, w sorted within
//! groups of equal axis lengths. Written as rows a.z <= b (z = (w, t)) plus a box |z_j| <= BOX, we
//! solve the dual   min b.y  s.t.  A^T y = 0, sum(y) = 1, y >= 0   by a dense revised simplex whose
//! basis is only (k+2) x (k+2). Its simplex multipliers (z, s) maximise the minimum slack s:
//! the system is feasible iff s > 0, and then z is a strictly separating witness.
//!
//! Usage:
//!   threshold-burnside <n> [--from m] [--id-known] [--threads T] [--out DIR]
//!   threshold-burnside grid <l1> <l2> ... [--threads T]

use std::collections::HashMap;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::{Duration, Instant};

// ---------------------------------------------------------------------------------------------
// Live progress (display only; exact counts are kept separately in Stats)
// ---------------------------------------------------------------------------------------------

static P_CANON: AtomicU64 = AtomicU64::new(0);
static P_LPS: AtomicU64 = AtomicU64::new(0);
static P_SUB_DONE: AtomicUsize = AtomicUsize::new(0);
static P_SUB_TOTAL: AtomicUsize = AtomicUsize::new(0);
static P_CLASS_IDX: AtomicUsize = AtomicUsize::new(0);
static P_CLASS_TOTAL: AtomicUsize = AtomicUsize::new(0);
static P_STOP: AtomicBool = AtomicBool::new(false);
static P_DONATED: AtomicU64 = AtomicU64::new(0);
static P_RUNNING: AtomicU64 = AtomicU64::new(0);
static P_IDLEG: AtomicU64 = AtomicU64::new(0);
static P_DCHECK: AtomicU64 = AtomicU64::new(0);
static P_DPEND: AtomicU64 = AtomicU64::new(0);
static P_DSHALLOW: AtomicU64 = AtomicU64::new(0);
static P_LABEL: Mutex<String> = Mutex::new(String::new());
static P_PHASE: Mutex<&'static str> = Mutex::new("starting");
static P_CLASS_START: Mutex<Option<Instant>> = Mutex::new(None);

fn progress_begin_class(label: String, idx: usize, total: usize) {
    *P_LABEL.lock().unwrap() = label;
    *P_PHASE.lock().unwrap() = "splitting search tree";
    *P_CLASS_START.lock().unwrap() = Some(Instant::now());
    P_CLASS_IDX.store(idx, Ordering::Relaxed);
    P_CLASS_TOTAL.store(total, Ordering::Relaxed);
    P_CANON.store(0, Ordering::Relaxed);
    P_LPS.store(0, Ordering::Relaxed);
    P_SUB_DONE.store(0, Ordering::Relaxed);
    P_SUB_TOTAL.store(0, Ordering::Relaxed);
}

fn sci(x: f64) -> String {
    if x < 1e5 { format!("{:.0}", x) } else { format!("{:.3e}", x) }
}

/// Prints one status line to stderr every `secs` seconds until P_STOP is set.
fn spawn_heartbeat(secs: u64) -> std::thread::JoinHandle<()> {
    std::thread::spawn(move || {
        let t0 = Instant::now();
        let (mut last_c, mut last_t) = (0u64, Instant::now());
        let mut next = Duration::from_secs(secs);
        while !P_STOP.load(Ordering::Relaxed) {
            std::thread::sleep(Duration::from_millis(200));
            if t0.elapsed() < next {
                continue;
            }
            next += Duration::from_secs(secs);
            let c = P_CANON.load(Ordering::Relaxed);
            if c < last_c {
                last_c = 0; // a new class started
            }
            let rate = (c - last_c) as f64 / last_t.elapsed().as_secs_f64().max(1e-9);
            last_c = c;
            last_t = Instant::now();
            let (done, tot) = (P_SUB_DONE.load(Ordering::Relaxed), P_SUB_TOTAL.load(Ordering::Relaxed));
            let cel = P_CLASS_START.lock().unwrap().map(|s| s.elapsed().as_secs_f64()).unwrap_or(0.0);
            let eta = if done > 0 && tot > done { format!("~{:.0}s", cel / done as f64 * (tot - done) as f64) } else { "?".into() };
            let sub = if tot > 0 { format!("{done}/{tot} ({:.1}%)", 100.0 * done as f64 / tot as f64) } else { "-".into() };
            eprintln!(
                "[{:>7.0}s] class {}/{} {} [{}]  subtrees {}  canonical {} ({}/s)  LPs {}  donated {}  class {:.0}s  ETA {}",
                t0.elapsed().as_secs_f64(),
                P_CLASS_IDX.load(Ordering::Relaxed),
                P_CLASS_TOTAL.load(Ordering::Relaxed),
                P_LABEL.lock().unwrap(),
                P_PHASE.lock().unwrap(),
                sub,
                sci(c as f64),
                sci(rate),
                sci(P_LPS.load(Ordering::Relaxed) as f64),
                format!("{} [running {} idle@lastpick {} checks {} pend/check {:.2} tooDeep {}]", sci(P_DONATED.load(Ordering::Relaxed) as f64), P_RUNNING.load(Ordering::Relaxed), P_IDLEG.load(Ordering::Relaxed), P_DCHECK.load(Ordering::Relaxed), P_DPEND.load(Ordering::Relaxed) as f64 / P_DCHECK.load(Ordering::Relaxed).max(1) as f64, P_DSHALLOW.load(Ordering::Relaxed)),
                cel,
                eta
            );
        }
    })
}

const MAXD: usize = 16; // max number of LP variables (k axes + threshold)
type Vd = [f64; MAXD];

const TOL: f64 = 1e-9;
const BOX: f64 = 1e4;
const FEAS: f64 = 1e-7; // max-min-slack must exceed this to count as feasible
const DECIDE: f64 = 1e-7; // |w.x - t| below this: witness is ambiguous, solve both branches
const MIN_TASKS: usize = 4096; // split the search tree until at least this many subtrees
const DONATE_MIN_REMAINING: usize = 6; // donate a branch to an idle worker only if >= this many points remain

/// Literature values of A000609(n) (all threshold functions of n variables) for the identity class.
fn a000609(n: usize) -> Option<u128> {
    [2u128, 4, 14, 104, 1882, 94572, 15028134, 8378070864, 17561539552946, 144130531453121108]
        .get(n)
        .copied()
}

/// A000617(n): NP-classes of threshold functions (= canonical count of the identity grid).
fn a000617(n: usize) -> Option<u64> {
    [2u64, 3, 5, 10, 27, 119, 1113, 29375, 2730166, 989913346].get(n).copied()
}

/// A002078(n): positive threshold functions (= positive count of the identity grid).
fn a002078(n: usize) -> Option<u128> {
    [2u128, 3, 6, 20, 150, 3287, 244158, 66291591, 68863243522].get(n).copied()
}

/// A001529(n): NPN-classes of threshold functions = (A000617(n) + #self-dual NP classes) / 2.
fn a001529(n: usize) -> Option<u64> {
    [1u64, 2, 3, 6, 15, 63, 567, 14755, 1366318].get(n).copied()
}

/// A001532(n): NP-classes of self-dual threshold functions (majority games) = self-dual canonical count.
/// (OEIS offset 1; n = 0 has none.)
fn a001532(n: usize) -> Option<u64> {
    [0u64, 1, 1, 2, 3, 7, 21, 135, 2470, 175428, 52980624].get(n).copied()
}

/// Report the identity-grid consistency checks and the NPN count derived from self-dual classes.
fn identity_report(n: usize, st: &Stats) {
    let ok = |b: bool| if b { "OK" } else { "MISMATCH" };
    if let Some(v) = a001532(n) {
        println!("    self-dual check vs A001532({n}) = {v}: {}", ok(v == st.selfdual));
    }
    if let Some(v) = a000609(n) {
        println!("    identity check vs A000609({n}) = {v}: {}", ok(v == st.total));
    }
    if let Some(v) = a000617(n) {
        println!("    canonical check vs A000617({n}) = {v}: {}", ok(v == st.canonical));
    }
    if let Some(v) = a002078(n) {
        println!("    positive check vs A002078({n}) = {v}: {}", ok(v == st.positive));
    }
    let twice = st.canonical + st.selfdual;
    let npn = twice / 2;
    let known = a001529(n).map(|v| format!("A001529({n}) = {v}: {}", ok(v == npn))).unwrap_or_else(|| "A001529 not known here: NEW VALUE".into());
    println!(
        "    NPN classes = (canonical {} + self-dual {}) / 2 = {}{}   [{}]",
        st.canonical,
        st.selfdual,
        npn,
        if twice % 2 == 1 { " (ODD SUM - ERROR)" } else { "" },
        known
    );
}

// ---------------------------------------------------------------------------------------------
// Grid
// ---------------------------------------------------------------------------------------------

struct Grid {
    l: Vec<usize>,
    k: usize,
    groups: Vec<Vec<usize>>,
    pf: Vec<Vd>,
    n: usize,
    preds: Vec<Vec<u32>>,
    steps: Vec<(Vec<u32>, Vec<u32>)>,
    swaps: Vec<Vec<Vec<u32>>>, // [group][adjacent pair] -> point permutation
    static_rows: Vec<(Vd, f64)>,
    h: u128,
    /// ps[i][a] = sum of coordinates of point i over axes b <= a in a's group. For weights that are
    /// >= 0 and non-increasing within groups: w.q <= w.p for all such w  iff  ps[q] <= ps[p] entrywise.
    ps: Vec<[i16; MAXD]>,
    comp: Vec<u32>, // index of the complementary point L - p
}

impl Grid {
    fn new(lengths: &[usize]) -> Grid {
        let mut l = lengths.to_vec();
        l.sort_unstable_by(|a, b| b.cmp(a));
        let k = l.len();
        assert!(k + 1 <= MAXD, "too many axes");
        let mut groups: Vec<Vec<usize>> = vec![];
        for a in 0..k {
            if a > 0 && l[a] == l[a - 1] {
                groups.last_mut().unwrap().push(a);
            } else {
                groups.push(vec![a]);
            }
        }
        let radix: Vec<usize> = l.iter().map(|x| x + 1).collect();
        let mut stride = vec![1usize; k];
        for a in (0..k.saturating_sub(1)).rev() {
            stride[a] = stride[a + 1] * radix[a + 1];
        }
        let total: usize = radix.iter().product();
        let mut pts: Vec<Vec<i32>> = (0..total)
            .map(|code| (0..k).map(|a| ((code / stride[a]) % radix[a]) as i32).collect())
            .collect();
        pts.sort_by(|p, q| {
            let sp: i32 = p.iter().sum();
            let sq: i32 = q.iter().sum();
            sp.cmp(&sq).then_with(|| p.cmp(q))
        });
        let enc = |p: &[i32]| -> usize { p.iter().zip(&stride).map(|(&x, &s)| x as usize * s).sum() };
        let mut idx = vec![0u32; total];
        for (i, p) in pts.iter().enumerate() {
            idx[enc(p)] = i as u32;
        }
        let mut preds = Vec::with_capacity(total);
        for p in &pts {
            let mut pr: Vec<u32> = vec![];
            for a in 0..k {
                if p[a] > 0 {
                    let mut q = p.clone();
                    q[a] -= 1;
                    pr.push(idx[enc(&q)]);
                }
            }
            for g in &groups {
                for x in 0..g.len() {
                    for &j in &g[x + 1..] {
                        let i = g[x];
                        if p[i] > 0 && (p[j] as usize) < l[j] {
                            let mut q = p.clone();
                            q[i] -= 1;
                            q[j] += 1;
                            pr.push(idx[enc(&q)]);
                        }
                    }
                }
            }
            pr.sort_unstable();
            pr.dedup();
            preds.push(pr);
        }
        let mut steps = vec![];
        for a in 0..k {
            let (mut fr, mut to) = (vec![], vec![]);
            for (i, p) in pts.iter().enumerate() {
                if (p[a] as usize) < l[a] {
                    let mut q = p.clone();
                    q[a] += 1;
                    fr.push(i as u32);
                    to.push(idx[enc(&q)]);
                }
            }
            steps.push((fr, to));
        }
        let mut swaps = vec![];
        for g in &groups {
            let mut gs = vec![];
            for w in g.windows(2) {
                let (a, b) = (w[0], w[1]);
                gs.push(
                    pts.iter()
                        .map(|p| {
                            let mut q = p.clone();
                            q.swap(a, b);
                            idx[enc(&q)]
                        })
                        .collect(),
                );
            }
            swaps.push(gs);
        }
        let comp: Vec<u32> = pts
            .iter()
            .map(|p| {
                let q: Vec<i32> = p.iter().zip(&l).map(|(&x, &len)| len as i32 - x).collect();
                idx[enc(&q)]
            })
            .collect();
        let mut pf = vec![[0.0; MAXD]; total];
        for (i, p) in pts.iter().enumerate() {
            for a in 0..k {
                pf[i][a] = p[a] as f64;
            }
        }
        // static LP rows over z = (w_0..w_{k-1}, t)
        let d = k + 1;
        let mut static_rows = vec![];
        for j in 0..d {
            let mut e = [0.0; MAXD];
            e[j] = 1.0;
            static_rows.push((e, BOX));
            let mut e = [0.0; MAXD];
            e[j] = -1.0;
            static_rows.push((e, BOX));
        }
        for j in 0..k {
            let mut e = [0.0; MAXD];
            e[j] = -1.0;
            static_rows.push((e, 0.0)); // w_j >= 0
        }
        for g in &groups {
            for w in g.windows(2) {
                let mut e = [0.0; MAXD];
                e[w[1]] = 1.0;
                e[w[0]] = -1.0;
                static_rows.push((e, 0.0)); // w_b <= w_a
            }
        }
        let h = groups.iter().map(|g| fact(g.len())).product();
        let mut ps = vec![[0i16; MAXD]; total];
        for (i, p) in pts.iter().enumerate() {
            for grp in &groups {
                let mut acc = 0i16;
                for &a in grp {
                    acc += p[a] as i16;
                    ps[i][a] = acc;
                }
            }
        }
        Grid { l, k, groups, pf, n: total, preds, steps, swaps, static_rows, h, ps, comp }
    }
}

fn fact(m: usize) -> u128 {
    (1..=m as u128).product()
}

// ---------------------------------------------------------------------------------------------
// Small dense LP: dual revised simplex (basis size d+1)
// ---------------------------------------------------------------------------------------------

/// Returns (s, z): max over z of min slack s, subject to rows a.z + s <= b. Feasible iff s > 0.
fn max_slack(rows: &[(Vd, f64)], d: usize) -> (f64, Vd) {
    let m = rows.len();
    let r = d + 1;
    let col = |j: usize, i: usize| -> f64 {
        if j < m {
            if i < d {
                rows[j].0[i]
            } else {
                1.0
            }
        } else if j - m == i {
            1.0
        } else {
            0.0
        }
    };
    // Fixed-size stack buffers: no heap allocation per LP (the Windows heap serialises threads).
    const R: usize = MAXD + 1;
    debug_assert!(r <= R);
    let mut basis_arr = [0usize; R];
    for i in 0..r {
        basis_arr[i] = m + i;
    }
    let basis = &mut basis_arr[..r];
    let mut binv = [0.0f64; R * R];
    let mut bmat = [0.0f64; R * R];
    let mut xb = [0.0f64; R];
    let mut pi = [0.0f64; R];
    let mut u = [0.0f64; R];
    for phase in 1..=2 {
        let cost = |j: usize| -> f64 {
            if phase == 1 {
                if j >= m { 1.0 } else { 0.0 }
            } else if j < m {
                rows[j].1
            } else {
                0.0
            }
        };
        let mut iter = 0usize;
        loop {
            iter += 1;
            assert!(iter < 20000, "simplex iteration limit");
            invert(&basis, &col, r, &mut bmat, &mut binv);
            for i in 0..r {
                xb[i] = binv[i * r + d]; // h = e_d
            }
            pi.iter_mut().for_each(|x| *x = 0.0);
            for i in 0..r {
                let cb = cost(basis[i]);
                if cb != 0.0 {
                    for t in 0..r {
                        pi[t] += cb * binv[i * r + t];
                    }
                }
            }
            let bland = iter > 60;
            let ncols = if phase == 1 { m + r } else { m };
            let mut enter = usize::MAX;
            let mut best = -TOL;
            for j in 0..ncols {
                if basis.contains(&j) {
                    continue;
                }
                let mut rj = cost(j);
                for i in 0..r {
                    rj -= pi[i] * col(j, i);
                }
                if rj < -TOL {
                    if bland {
                        enter = j;
                        break;
                    }
                    if rj < best {
                        best = rj;
                        enter = j;
                    }
                }
            }
            if enter == usize::MAX {
                if phase == 1 {
                    let infeas: f64 = (0..r).filter(|&i| basis[i] >= m).map(|i| xb[i]).sum();
                    assert!(infeas < 1e-7, "phase 1 failed: {infeas}");
                    // drive remaining artificials out of the basis where possible
                    for i in 0..r {
                        if basis[i] < m {
                            continue;
                        }
                        for j in 0..m {
                            if basis.contains(&j) {
                                continue;
                            }
                            let v: f64 = (0..r).map(|t| binv[i * r + t] * col(j, t)).sum();
                            if v.abs() > 1e-7 {
                                basis[i] = j;
                                invert(&basis, &col, r, &mut bmat, &mut binv);
                                break;
                            }
                        }
                    }
                    break;
                } else {
                    let mut z = [0.0; MAXD];
                    z[..d].copy_from_slice(&pi[..d]);
                    return (pi[d], z);
                }
            }
            for i in 0..r {
                u[i] = (0..r).map(|t| binv[i * r + t] * col(enter, t)).sum();
            }
            let mut leave = usize::MAX;
            let mut best_ratio = f64::INFINITY;
            for i in 0..r {
                if u[i] > TOL {
                    let ratio = xb[i].max(0.0) / u[i];
                    if ratio < best_ratio - 1e-12
                        || (ratio <= best_ratio + 1e-12 && leave != usize::MAX && basis[i] < basis[leave])
                    {
                        best_ratio = ratio;
                        leave = i;
                    }
                }
            }
            assert!(leave != usize::MAX, "LP unbounded (impossible: y lies in a simplex)");
            basis[leave] = enter;
        }
    }
    unreachable!()
}

/// binv = inverse of the matrix whose columns are the basis columns (Gauss-Jordan, partial pivoting).
fn invert(basis: &[usize], col: &dyn Fn(usize, usize) -> f64, r: usize, a: &mut [f64], inv: &mut [f64]) {
    for i in 0..r {
        for t in 0..r {
            a[i * r + t] = col(basis[t], i);
            inv[i * r + t] = if i == t { 1.0 } else { 0.0 };
        }
    }
    for c in 0..r {
        let mut p = c;
        for i in c + 1..r {
            if a[i * r + c].abs() > a[p * r + c].abs() {
                p = i;
            }
        }
        assert!(a[p * r + c].abs() > 1e-12, "singular basis");
        if p != c {
            for t in 0..r {
                a.swap(p * r + t, c * r + t);
                inv.swap(p * r + t, c * r + t);
            }
        }
        let piv = a[c * r + c];
        for t in 0..r {
            a[c * r + t] /= piv;
            inv[c * r + t] /= piv;
        }
        for i in 0..r {
            if i != c {
                let f = a[i * r + c];
                if f != 0.0 {
                    for t in 0..r {
                        a[i * r + t] -= f * a[c * r + t];
                        inv[i * r + t] -= f * inv[c * r + t];
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Search
// ---------------------------------------------------------------------------------------------

#[derive(Clone, Copy, Default, Debug)]
struct Stats {
    total: u128,
    positive: u128,
    canonical: u64,
    lps: u64,
    bad: u64,
    selfdual: u64, // canonical functions with f(L - x) = not f(x): self-dual (NPN-self-complementary on the cube)
}

impl Stats {
    fn add(mut self, o: Stats) -> Stats {
        self.total += o.total;
        self.positive += o.positive;
        self.canonical += o.canonical;
        self.lps += o.lps;
        self.bad += o.bad;
        self.selfdual += o.selfdual;
        self
    }
}

/// A branch point on the explicit DFS stack.
struct Frame {
    i: usize,
    #[allow(dead_code)]
    cur: i8,
    pending: Option<(i8, Vd)>, // unexplored alternative branch and its witness
    n_ones: usize,             // lengths of ones/zeros before this point was labelled
    n_zeros: usize,
}

#[derive(Clone)]
struct Task {
    i: usize,
    lab: Vec<i8>,
    ones: Vec<u32>,
    zeros: Vec<u32>,
    z: Vd,
}

struct Search<'a> {
    g: &'a Grid,
    lab: Vec<i8>,
    ones: Vec<u32>,
    zeros: Vec<u32>,
    rows: Vec<(Vd, f64)>,
    split: Option<usize>,
    tasks: Vec<Task>,
    st: Stats,
    pend_canon: u64,
    pend_lps: u64,
    ctx: Option<&'a Shared<'a>>,
    root: usize,
    zbuf: Vec<u32>,
    fbuf: Vec<bool>,
}

impl Drop for Search<'_> {
    fn drop(&mut self) {
        self.flush();
    }
}

impl<'a> Search<'a> {
    fn new(g: &'a Grid) -> Self {
        Search {
            g,
            lab: vec![-1; g.n],
            ones: vec![],
            zeros: vec![],
            rows: Vec::with_capacity(g.n + g.static_rows.len()),
            split: None,
            tasks: vec![],
            st: Stats::default(),
            pend_canon: 0,
            pend_lps: 0,
            ctx: None,
            root: 0,
            zbuf: Vec::new(),
            fbuf: Vec::new(),
        }
    }

    fn flush(&mut self) {
        P_CANON.fetch_add(self.pend_canon, Ordering::Relaxed);
        P_LPS.fetch_add(self.pend_lps, Ordering::Relaxed);
        self.pend_canon = 0;
        self.pend_lps = 0;
    }

    fn lp(&mut self) -> Option<Vd> {
        let g = self.g;
        let k = g.k;
        self.rows.clear();
        self.rows.extend_from_slice(&g.static_rows);
        for &i in &self.ones {
            let mut a = [0.0; MAXD];
            for j in 0..k {
                a[j] = -g.pf[i as usize][j];
            }
            a[k] = 1.0;
            self.rows.push((a, 0.0));
        }
        // Only maximal zeros matter: a zero dominated by a later zero is implied by it.
        // (Zeros are pushed in linear-extension order, so only later ones can dominate earlier ones.)
        self.zbuf.clear();
        for (x, &q) in self.zeros.iter().enumerate() {
            let pq = &g.ps[q as usize];
            let dominated = self.zeros[x + 1..].iter().any(|&p| {
                let pp = &g.ps[p as usize];
                (0..k).all(|a| pq[a] <= pp[a])
            });
            if !dominated {
                self.zbuf.push(q);
            }
        }
        for &i in &self.zbuf {
            let mut a = [0.0; MAXD];
            for j in 0..k {
                a[j] = g.pf[i as usize][j];
            }
            a[k] = -1.0;
            self.rows.push((a, -1.0));
        }
        self.st.lps += 1;
        self.pend_lps += 1;
        if self.pend_lps >= 1024 {
            self.flush();
        }
        let (s, z) = max_slack(&self.rows, k + 1);
        if s > FEAS { Some(z) } else { None }
    }

    fn value(&self, i: usize, z: &Vd) -> f64 {
        let k = self.g.k;
        let p = &self.g.pf[i];
        let mut v = -z[k];
        for j in 0..k {
            v += p[j] * z[j];
        }
        v
    }

    fn start(&mut self, split: Option<usize>) -> Vec<Task> {
        self.split = split;
        let z = self.lp().expect("empty system infeasible?");
        self.dfs(0, z, 0);
        std::mem::take(&mut self.tasks)
    }

    /// Continue from a task, keeping the current `split` setting (used for iterative deepening).
    fn resume_split(&mut self, t: &Task) {
        self.lab.copy_from_slice(&t.lab);
        self.ones = t.ones.clone();
        self.zeros = t.zeros.clone();
        self.dfs(t.i, t.z, 0);
    }

    fn resume(&mut self, t: &Task) {
        self.split = None;
        self.lab.copy_from_slice(&t.lab);
        self.ones = t.ones.clone();
        self.zeros = t.zeros.clone();
        self.run_iter(t.i, t.z);
    }

    /// Witness for each branch at point i: the free branch costs nothing, the other needs an LP.
    fn branch_witnesses(&mut self, i: usize, z: &Vd) -> [(i8, Option<Vd>); 2] {
        let v = self.value(i, z);
        let free: Option<i8> = if v >= DECIDE { Some(1) } else if v <= -DECIDE { Some(0) } else { None };
        let order: [i8; 2] = if free == Some(0) { [0, 1] } else { [1, 0] };
        let mut feas: [(i8, Option<Vd>); 2] = [(order[0], None), (order[1], None)];
        for slot in feas.iter_mut() {
            let br = slot.0;
            slot.1 = if Some(br) == free {
                Some(*z)
            } else {
                if br == 1 { self.ones.push(i as u32) } else { self.zeros.push(i as u32) }
                let r = self.lp();
                if br == 1 { self.ones.pop(); } else { self.zeros.pop(); }
                r
            };
        }
        feas
    }

    /// Iterative DFS with an explicit stack of branch points. When a worker is idle, the
    /// *shallowest* unexplored alternative (the biggest remaining piece of work) is donated.
    fn run_iter(&mut self, i0: usize, z0: Vd) {
        let g = self.g;
        let mut stack: Vec<Frame> = Vec::new();
        let (mut i, mut z) = (i0, z0);
        loop {
            // descend to a leaf
            loop {
                if let Some(sh) = self.ctx {
                    if sh.idle.load(Ordering::Relaxed) > 0 {
                        self.donate_shallowest(&mut stack, sh);
                    }
                }
                while i < g.n && g.preds[i].iter().any(|&q| self.lab[q as usize] == 1) {
                    self.lab[i] = 1;
                    i += 1;
                }
                if i == g.n {
                    self.leaf(&z);
                    break;
                }
                let feas = self.branch_witnesses(i, &z);
                let mut it = feas.into_iter().filter_map(|(b, w)| w.map(|w| (b, w)));
                let (b0, z0) = it.next().expect("no feasible branch below a feasible node");
                let pending = it.next();
                stack.push(Frame { i, cur: b0, pending, n_ones: self.ones.len(), n_zeros: self.zeros.len() });
                self.lab[i] = b0;
                if b0 == 1 { self.ones.push(i as u32) } else { self.zeros.push(i as u32) }
                z = z0;
                i += 1;
            }
            // backtrack to the deepest frame with an unexplored alternative
            loop {
                let Some(f) = stack.pop() else { return };
                for x in &mut self.lab[f.i..] {
                    *x = -1;
                }
                self.ones.truncate(f.n_ones);
                self.zeros.truncate(f.n_zeros);
                if let Some((b, zz)) = f.pending {
                    self.lab[f.i] = b;
                    if b == 1 { self.ones.push(f.i as u32) } else { self.zeros.push(f.i as u32) }
                    stack.push(Frame { i: f.i, cur: b, pending: None, n_ones: f.n_ones, n_zeros: f.n_zeros });
                    i = f.i + 1;
                    z = zz;
                    break;
                }
            }
        }
    }

    fn donate_shallowest(&mut self, stack: &mut [Frame], sh: &Shared) {
        let g = self.g;
        if std::env::var_os("TB_DEBUG").is_some() {
            P_DCHECK.fetch_add(1, Ordering::Relaxed);
            let np = stack.iter().filter(|f| f.pending.is_some()).count() as u64;
            P_DPEND.fetch_add(np, Ordering::Relaxed);
            if stack.iter().any(|f| f.pending.is_some() && g.n - f.i < DONATE_MIN_REMAINING) { P_DSHALLOW.fetch_add(1, Ordering::Relaxed); }
        }
        for f in stack.iter_mut() {
            if f.pending.is_some() && g.n - f.i >= DONATE_MIN_REMAINING {
                let (b, zz) = f.pending.take().unwrap();
                // state as it was when this frame was created: lab[..f.i] is unchanged since then
                let mut lab = vec![-1i8; g.n];
                lab[..f.i].copy_from_slice(&self.lab[..f.i]);
                lab[f.i] = b;
                let mut ones = self.ones[..f.n_ones].to_vec();
                let mut zeros = self.zeros[..f.n_zeros].to_vec();
                if b == 1 { ones.push(f.i as u32) } else { zeros.push(f.i as u32) }
                sh.donate(self.root, Task { i: f.i + 1, lab, ones, zeros, z: zz });
                return;
            }
        }
    }

    fn dfs(&mut self, mut i: usize, z: Vd, depth: usize) {
        let g = self.g;
        if self.split == Some(depth) {
            self.tasks.push(Task { i, lab: self.lab.clone(), ones: self.ones.clone(), zeros: self.zeros.clone(), z });
            return;
        }
        while i < g.n && g.preds[i].iter().any(|&q| self.lab[q as usize] == 1) {
            self.lab[i] = 1;
            i += 1;
        }
        if i == g.n {
            self.leaf(&z);
            return;
        }
        let v = self.value(i, &z);
        let free: Option<i8> = if v >= DECIDE { Some(1) } else if v <= -DECIDE { Some(0) } else { None };
        let order: [i8; 2] = match free {
            Some(1) => [1, 0],
            Some(_) => [0, 1],
            None => [1, 0],
        };
        // Witnesses for the feasible branches (the free one costs nothing; the other needs an LP).
        let mut feas: [(i8, Option<Vd>); 2] = [(order[0], None), (order[1], None)];
        for slot in feas.iter_mut() {
            let br = slot.0;
            slot.1 = if Some(br) == free {
                Some(z)
            } else {
                if br == 1 { self.ones.push(i as u32) } else { self.zeros.push(i as u32) }
                let r = self.lp();
                if br == 1 { self.ones.pop(); } else { self.zeros.pop(); }
                r
            };
        }
        for (br, zz) in feas {
            if let Some(zz) = zz {
                self.branch(i, br, zz, depth);
            }
        }
    }

    /// Label point i with `br` and search the subtree below it; restores the state afterwards.
    fn branch(&mut self, i: usize, br: i8, z: Vd, depth: usize) {
        if br == 1 { self.ones.push(i as u32) } else { self.zeros.push(i as u32) }
        self.lab[i] = br;
        self.dfs(i + 1, z, depth + 1);
        for x in &mut self.lab[i..] {
            *x = -1;
        }
        if br == 1 { self.ones.pop(); } else { self.zeros.pop(); }
    }

    fn leaf(&mut self, z: &Vd) {
        let g = self.g;
        let mut f = std::mem::take(&mut self.fbuf);
        f.clear();
        f.extend(self.lab.iter().map(|&x| x == 1));
        for i in 0..g.n {
            let v = self.value(i, z);
            if (f[i] && v < -TOL) || (!f[i] && v >= -TOL) {
                self.st.bad += 1;
                break;
            }
        }
        let mut ess = 0u32;
        for a in 0..g.k {
            let (fr, to) = &g.steps[a];
            if fr.iter().zip(to).any(|(&x, &y)| f[x as usize] != f[y as usize]) {
                ess += 1;
            }
        }
        let mut stab: u128 = 1;
        for (gi, grp) in g.groups.iter().enumerate() {
            let mut run = 1usize;
            for pair in 0..grp.len().saturating_sub(1) {
                let perm = &g.swaps[gi][pair];
                if (0..g.n).all(|x| f[x] == f[perm[x] as usize]) {
                    run += 1;
                } else {
                    stab *= fact(run);
                    run = 1;
                }
            }
            stab *= fact(run);
        }
        let orbit = g.h / stab;
        if (0..g.n).all(|x| f[x] != f[g.comp[x] as usize]) {
            self.st.selfdual += 1;
        }
        self.st.canonical += 1;
        self.pend_canon += 1;
        if self.pend_canon >= 1024 {
            self.flush();
        }
        self.st.positive += orbit;
        self.st.total += orbit << ess;
        self.fbuf = f;
    }
}

// ---------------------------------------------------------------------------------------------
// Grid counting with task splitting, parallelism and checkpoints
// ---------------------------------------------------------------------------------------------

fn lam_name(l: &[usize]) -> String {
    l.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("-")
}

fn count_grid(lengths: &[usize], ckpt_dir: Option<&Path>, verbose: bool) -> (Stats, f64, usize) {
    let t0 = Instant::now();
    let g = Grid::new(lengths);
    // Split the tree by iterative deepening: each round expands every current subtree by 2 more
    // branch levels (leaves met on the way are counted once, in `pre`). Deterministic and
    // independent of the thread count, so checkpoints stay valid across runs.
    let mut depth = 4;
    let mut s0 = Search::new(&g);
    let mut tasks = s0.start(Some(depth));
    let mut pre = s0.st;
    drop(s0);
    while !tasks.is_empty() && tasks.len() < MIN_TASKS && depth < 64 {
        let mut next = vec![];
        for t in &tasks {
            let mut s = Search::new(&g);
            s.split = Some(2);
            s.resume_split(t);
            pre = pre.add(s.st);
            next.append(&mut s.tasks);
        }
        tasks = next;
        depth += 2;
    }
    let ntasks = tasks.len();
    P_SUB_TOTAL.store(ntasks, Ordering::Relaxed);
    *P_PHASE.lock().unwrap() = "searching";
    // checkpoint: one line per finished task
    let mut done: HashMap<usize, Stats> = HashMap::new();
    let ck_path: Option<PathBuf> =
        ckpt_dir.map(|d| d.join(format!("grid_{}_d{}_t{}.tasks", lam_name(&g.l), depth, ntasks)));
    if let Some(p) = &ck_path {
        if let Ok(f) = fs::File::open(p) {
            for line in BufReader::new(f).lines().map_while(Result::ok) {
                let v: Vec<&str> = line.split_whitespace().collect();
                if v.len() == 6 || v.len() == 7 {
                    done.insert(
                        v[0].parse().unwrap(),
                        Stats {
                            total: v[1].parse().unwrap(),
                            positive: v[2].parse().unwrap(),
                            canonical: v[3].parse().unwrap(),
                            lps: v[4].parse().unwrap(),
                            bad: v[5].parse().unwrap(),
                            selfdual: v.get(6).map(|x| x.parse().unwrap()).unwrap_or(0),
                        },
                    );
                }
            }
        }
    }
    let resumed = done.len();
    let writer = ck_path.as_ref().map(|p| Mutex::new(OpenOptions::new().create(true).append(true).open(p).unwrap()));
    P_SUB_DONE.store(resumed, Ordering::Relaxed);
    if verbose && resumed > 0 {
        eprintln!("    [{}] resuming: {}/{} subtrees already done (checkpoint)", lam_name(&g.l), resumed, ntasks);
        P_CANON.fetch_add(done.values().map(|s| s.canonical).sum::<u64>(), Ordering::Relaxed);
    }
    let todo: Vec<usize> = (0..ntasks).filter(|i| !done.contains_key(i)).collect();
    // Shared work queue. Each checkpoint task is the root of a "family": work donated from it
    // stays in the family, and the root is checkpointed once the whole family has finished.
    let sh = Shared {
        g: &g,
        queue: Mutex::new(todo.iter().rev().map(|&ti| (ti, tasks[ti].clone())).collect()),
        cv: Condvar::new(),
        idle: AtomicUsize::new(0),
        active: AtomicUsize::new(todo.len()),
        fam: (0..ntasks).map(|ti| AtomicUsize::new(if done.contains_key(&ti) { 0 } else { 1 })).collect(),
        fam_stats: (0..ntasks).map(|_| Mutex::new(Stats::default())).collect(),
        writer,
        new: Mutex::new(Stats::default()),
    };
    drop(tasks);
    let nthreads = num_threads();
    std::thread::scope(|scope| {
        for _ in 0..nthreads {
            scope.spawn(|| sh.worker());
        }
    });
    *P_PHASE.lock().unwrap() = "finished";
    let mut st = pre.add(*sh.new.lock().unwrap());
    for s in done.values() {
        st = st.add(*s);
    }
    (st, t0.elapsed().as_secs_f64(), ntasks)
}

static THREADS: AtomicUsize = AtomicUsize::new(0);

fn num_threads() -> usize {
    match THREADS.load(Ordering::Relaxed) {
        0 => std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        t => t,
    }
}

struct Shared<'a> {
    g: &'a Grid,
    queue: Mutex<Vec<(usize, Task)>>, // (root task index, task); LIFO so donated work runs first
    cv: Condvar,
    idle: AtomicUsize,        // workers waiting for work
    active: AtomicUsize,      // queued or running items
    fam: Vec<AtomicUsize>,    // live items per root
    fam_stats: Vec<Mutex<Stats>>,
    writer: Option<Mutex<fs::File>>,
    new: Mutex<Stats>,        // total over roots finished in this run
}

impl<'a> Shared<'a> {
    fn donate(&self, root: usize, t: Task) {
        self.fam[root].fetch_add(1, Ordering::SeqCst);
        let mut q = self.queue.lock().unwrap();
        self.active.fetch_add(1, Ordering::SeqCst);
        q.push((root, t));
        P_DONATED.fetch_add(1, Ordering::Relaxed);
        drop(q);
        self.cv.notify_one();
    }

    fn worker(&self) {
        loop {
            let item = {
                let mut q = self.queue.lock().unwrap();
                loop {
                    if let Some(it) = q.pop() {
                        break Some(it);
                    }
                    if self.active.load(Ordering::SeqCst) == 0 {
                        break None;
                    }
                    self.idle.fetch_add(1, Ordering::SeqCst);
                    q = self.cv.wait(q).unwrap();
                    self.idle.fetch_sub(1, Ordering::SeqCst);
                }
            };
            let Some((root, task)) = item else {
                self.cv.notify_all();
                return;
            };
            let mut s = Search::new(self.g);
            s.ctx = Some(self);
            s.root = root;
            P_RUNNING.fetch_add(1, Ordering::Relaxed);
            P_IDLEG.store(self.idle.load(Ordering::Relaxed) as u64, Ordering::Relaxed);
            let t_item = Instant::now();
            let depth_i = task.i;
            s.resume(&task);
            let st = std::mem::take(&mut s.st);
            P_RUNNING.fetch_sub(1, Ordering::Relaxed);
            if std::env::var("TB_DEBUG").is_ok() && t_item.elapsed().as_secs_f64() > 2.0 {
                eprintln!("      slow item: root {root} start i={depth_i}/{} canonical={} lps={} {:.1}s",
                    self.g.n, st.canonical, st.lps, t_item.elapsed().as_secs_f64());
            }
            drop(s);
            {
                let mut fs_ = self.fam_stats[root].lock().unwrap();
                *fs_ = fs_.add(st);
            }
            if self.fam[root].fetch_sub(1, Ordering::SeqCst) == 1 {
                let fst = *self.fam_stats[root].lock().unwrap();
                if let Some(w) = &self.writer {
                    writeln!(w.lock().unwrap(), "{} {} {} {} {} {} {}", root, fst.total, fst.positive, fst.canonical, fst.lps, fst.bad, fst.selfdual)
                        .unwrap();
                }
                {
                    let mut nw = self.new.lock().unwrap();
                    *nw = nw.add(fst);
                }
                P_SUB_DONE.fetch_add(1, Ordering::Relaxed);
            }
            let q = self.queue.lock().unwrap();
            if self.active.fetch_sub(1, Ordering::SeqCst) == 1 {
                drop(q);
                self.cv.notify_all();
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Burnside driver
// ---------------------------------------------------------------------------------------------

fn partitions(n: usize, maxp: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if n == 0 {
        out.push(cur.clone());
        return;
    }
    for p in (1..=n.min(maxp)).rev() {
        cur.push(p);
        partitions(n - p, p, cur, out);
        cur.pop();
    }
}

fn class_size(lam: &[usize]) -> u128 {
    let n: usize = lam.iter().sum();
    let mut cnt: HashMap<usize, usize> = HashMap::new();
    for &l in lam {
        *cnt.entry(l).or_default() += 1;
    }
    let mut d: u128 = 1;
    for (&l, &m) in &cnt {
        d *= (l as u128).pow(m as u32) * fact(m);
    }
    fact(n) / d
}

fn run_n(n: usize, id_known: bool, out: &Path) {
    let t0 = Instant::now();
    fs::create_dir_all(out).unwrap();
    let res_path = out.join(format!("n{n}.tsv"));
    // resume finished classes
    let mut have: HashMap<String, u128> = HashMap::new();
    if let Ok(f) = fs::File::open(&res_path) {
        for line in BufReader::new(f).lines().map_while(Result::ok) {
            let v: Vec<&str> = line.split('\t').collect();
            if v.len() >= 2 && v[0] != "lambda" {
                have.insert(v[0].to_string(), v[1].parse().unwrap());
            }
        }
    }
    let new_file = !res_path.exists();
    let mut rf = OpenOptions::new().create(true).append(true).open(&res_path).unwrap();
    if new_file {
        writeln!(rf, "lambda\tfix\tcanonical\tpositive\tlps\tbad\tsecs\tsubtrees\tselfdual").unwrap();
    }
    let mut lams = vec![];
    partitions(n, n, &mut vec![], &mut lams);
    lams.sort_by_key(|l| l.len()); // cheap grids first
    println!("n={n}:");
    let mut sum: u128 = 0;
    for (ci, lam) in lams.iter().enumerate() {
        let name = lam_name(lam);
        progress_begin_class(format!("n={n} {name}"), ci + 1, lams.len());
        let cs = class_size(lam);
        let is_id = lam.iter().all(|&x| x == 1);
        let fix = if is_id && id_known {
            let v = a000609(n).expect("A000609(n) unknown: identity term must be computed");
            println!("  {name}: class {cs}, Fix {v}  (A000609, literature)");
            v
        } else if let Some(&v) = have.get(&name) {
            println!("  {name}: class {cs}, Fix {v}  (from {})", res_path.display());
            v
        } else {
            eprintln!("  -> starting class {}/{}: {name} (class size {cs})", ci + 1, lams.len());
            let (st, secs, nt) = count_grid(lam, Some(out), true);
            println!(
                "  {name}: class {cs}, Fix {}  canonical={} positive={} selfdual={} lps={} bad={} {:.1}s subtrees={}",
                st.total, st.canonical, st.positive, st.selfdual, st.lps, st.bad, secs, nt
            );
            // selfdual is a trailing column (files created before it existed keep the old header)
            writeln!(
                rf,
                "{name}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                st.total, st.canonical, st.positive, st.lps, st.bad, secs, nt, st.selfdual
            )
            .unwrap();
            if is_id {
                identity_report(n, &st);
            }
            st.total
        };
        sum += cs * fix;
    }
    let nf = fact(n);
    println!("==> a({n}) = {}   (Burnside remainder {}, {:.1}s)", sum / nf, sum % nf, t0.elapsed().as_secs_f64());
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
    let hb = if every > 0 { Some(spawn_heartbeat(every)) } else { None };
    if args.first().map(String::as_str) == Some("grid") {
        let ls: Vec<usize> = args[1..].iter().take_while(|a| !a.starts_with("--")).map(|a| a.parse().unwrap()).collect();
        progress_begin_class(format!("grid {}", lam_name(&ls)), 1, 1);
        let (st, secs, nt) = count_grid(&ls, None, true);
        println!(
            "grid {:?}: total={} positive={} canonical={} selfdual={} lps={} bad={} {:.2}s subtrees={}",
            ls, st.total, st.positive, st.canonical, st.selfdual, st.lps, st.bad, secs, nt
        );
        if ls.iter().all(|&x| x == 1) {
            identity_report(ls.len(), &st);
        }
    } else {
        let n_max: usize = args.first().expect("usage: threshold-burnside <n> | grid <l...>").parse().unwrap();
        let n_min: usize = opt("--from").map(|s| s.parse().unwrap()).unwrap_or(n_max);
        let id_known = args.iter().any(|a| a == "--id-known");
        let out = PathBuf::from(opt("--out").unwrap_or_else(|| "results".into()));
        for n in n_min..=n_max {
            run_n(n, id_known, &out);
        }
    }
    P_STOP.store(true, Ordering::Relaxed);
    if let Some(h) = hb {
        h.join().ok();
    }
}
