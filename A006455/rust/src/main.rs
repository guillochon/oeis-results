//! A006455: naturally labeled posets on [N], via the split identity
//!
//!   a(n+k) = sum over unlabeled P on n points, S on k points of w(P) w(S) Hom(S, J(P)),
//!
//! where w(P) = e(P)/|Aut P| (number of natural labelings of P), J(P) is the lattice of order
//! ideals of P and Hom(S, J(P)) counts order-preserving maps S -> J(P) (= ideals of S^op x P).
//! See ../README.md for the proof. Posets come from nauty's genposetg (../gen_posets.sh).
//!
//! Usage:
//!   a006455 check [nmax]        per-poset data: sum w = A006455(n), sum n!/|Aut| = A001035(n)
//!   a006455 split n k [threads] [flags]
//!                               a(n+k) from posets on n (lattice side) and k (structure side);
//!                               flags (default ods): o = optimized processing order, d = cheaper
//!                               direction per pair, s = unordered pairs when n = k; - = none

use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const MAXN: usize = 12;

/// A006455(0..12), for checks.
const A006455: [u128; 13] = [
    1, 1, 2, 7, 40, 357, 4824, 96428, 2800472, 116473461, 6855780268, 565505147444, 64824245807684,
];
/// A001035(0..10): labeled posets.
const A001035: [u128; 11] =
    [1, 1, 3, 19, 219, 4231, 130023, 6129859, 431723379, 44511042511, 6611065248783];

/// A poset in natural labeling: down[j] is the strict down-set of j, a subset of {0..j-1}.
#[derive(Clone)]
struct Poset {
    n: usize,
    down: [u16; MAXN],
    /// Hasse diagram: cov[j] = elements covered by j.
    cov: [u16; MAXN],
}

/// Parse one digraph6 line from genposetg's `t` output (edges x->y have x<y; Hasse diagram).
fn parse_d6(line: &str) -> Poset {
    let b = line.trim().as_bytes();
    assert!(b[0] == b'&', "not digraph6: {line}");
    let n = (b[1] - 63) as usize;
    assert!(n <= MAXN);
    let bit = |p: usize| -> bool {
        let c = (b[2 + p / 6] - 63) as u32;
        c >> (5 - p % 6) & 1 == 1
    };
    let mut cov = [0u16; MAXN];
    for x in 0..n {
        for y in 0..n {
            if bit(x * n + y) {
                assert!(x < y, "edge {x}->{y} not topological");
                cov[y] |= 1 << x;
            }
        }
    }
    let mut down = [0u16; MAXN];
    for y in 0..n {
        let mut d = 0u16;
        for x in 0..y {
            if cov[y] >> x & 1 == 1 {
                d |= down[x] | 1 << x;
            }
        }
        down[y] = d;
    }
    Poset { n, down, cov }
}

fn load(n: usize) -> Vec<Poset> {
    let path = format!("{}/../data/posets{n}.d6", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e} (run gen_posets.sh)"));
    text.lines().filter(|l| !l.is_empty()).map(parse_d6).collect()
}

/// Order ideals of p as bitmasks (unsorted).
fn ideals(p: &Poset) -> Vec<u16> {
    let mut res = vec![0u16];
    for j in 0..p.n {
        let d = p.down[j];
        let len = res.len();
        for t in 0..len {
            let m = res[t];
            if m & d == d {
                res.push(m | 1 << j);
            }
        }
    }
    res
}

/// Number of linear extensions: DP over ideals (maximal element removed last).
fn linear_extensions(p: &Poset) -> u128 {
    let mut ids = ideals(p);
    ids.sort_by_key(|m| m.count_ones());
    let mut ways = vec![0u128; 1 << p.n];
    let mut is_ideal = vec![false; 1 << p.n];
    for &m in &ids {
        is_ideal[m as usize] = true;
    }
    ways[0] = 1;
    for &m in &ids[1..] {
        let mut s = 0;
        let mut r = m;
        while r != 0 {
            let x = r.trailing_zeros();
            r &= r - 1;
            let q = m ^ 1 << x;
            if is_ideal[q as usize] {
                s += ways[q as usize];
            }
        }
        ways[m as usize] = s;
    }
    ways[(1usize << p.n) - 1]
}

/// |Aut p| by backtracking over relation-preserving bijections.
fn automorphisms(p: &Poset) -> u128 {
    let n = p.n;
    let mut up = [0u16; MAXN];
    for y in 0..n {
        for x in 0..n {
            if p.down[y] >> x & 1 == 1 {
                up[x] |= 1 << y;
            }
        }
    }
    let sig: Vec<(u32, u32, u32, u32)> = (0..n)
        .map(|x| (p.down[x].count_ones(), up[x].count_ones(), p.cov[x].count_ones(), 0))
        .collect();
    fn rec(i: usize, n: usize, p: &Poset, sig: &[(u32, u32, u32, u32)], img: &mut [usize; MAXN], used: u16) -> u128 {
        if i == n {
            return 1;
        }
        let mut total = 0;
        for v in 0..n {
            if used >> v & 1 == 1 || sig[v] != sig[i] {
                continue;
            }
            let ok = (0..i).all(|j| {
                let a = img[j];
                (p.down[i] >> j & 1) == (p.down[v] >> a & 1) && (p.down[j] >> i & 1) == (p.down[a] >> v & 1)
            });
            if ok {
                img[i] = v;
                total += rec(i + 1, n, p, sig, img, used | 1 << v);
            }
        }
        total
    }
    let mut img = [0usize; MAXN];
    rec(0, n, p, &sig, &mut img, 0)
}

/// Natural labelings of p up to equality: e(p)/|Aut p| (Aut acts freely on linear extensions).
fn weight(p: &Poset) -> u128 {
    let e = linear_extensions(p);
    let a = automorphisms(p);
    assert!(e % a == 0);
    e / a
}

/// The same poset with element order[t] renamed t (order must be a linear extension).
fn relabel(p: &Poset, order: &[usize]) -> Poset {
    let n = p.n;
    let mut new_of = [0usize; MAXN];
    for (t, &x) in order.iter().enumerate() {
        new_of[x] = t;
    }
    let map = |m: u16| -> u16 { (0..n).filter(|&x| m >> x & 1 == 1).fold(0, |a, x| a | 1 << new_of[x]) };
    let mut q = Poset { n, down: [0; MAXN], cov: [0; MAXN] };
    for x in 0..n {
        q.down[new_of[x]] = map(p.down[x]);
        q.cov[new_of[x]] = map(p.cov[x]);
    }
    for t in 0..n {
        assert!(q.down[t] >> t == 0, "relabel: not a linear extension");
    }
    q
}

/// The dual poset, naturally labeled by reversing labels.
fn dual(p: &Poset) -> Poset {
    let n = p.n;
    let mut q = Poset { n, down: [0; MAXN], cov: [0; MAXN] };
    for y in 0..n {
        for x in 0..n {
            if p.down[y] >> x & 1 == 1 {
                q.down[n - 1 - x] |= 1 << (n - 1 - y);
            }
            if p.cov[y] >> x & 1 == 1 {
                q.cov[n - 1 - x] |= 1 << (n - 1 - y);
            }
        }
    }
    q
}

/// Number of DP classes once the elements of ideal `i` are processed: active elements (with a
/// coverer outside i) grouped by their set of coverers outside i.
fn classes_after(p: &Poset, i: u16) -> u32 {
    let mut keys: Vec<u16> = Vec::new();
    for x in 0..p.n {
        if i >> x & 1 == 0 {
            continue;
        }
        let out = (0..p.n).filter(|&t| i >> t & 1 == 0 && p.cov[t] >> x & 1 == 1).fold(0u16, |a, t| a | 1 << t);
        if out != 0 && !keys.contains(&out) {
            keys.push(out);
        }
    }
    keys.len() as u32
}

/// DP cost model for one step over a lattice of size l: scan the old states (each scanning about
/// l/2 upper bounds when j keeps a class), then clear the new state array.
fn step_cost(l: f64, c_old: u32, has_j: bool, c_new: u32) -> f64 {
    l.powi(c_old as i32) * if has_j { 0.5 * l } else { 1.0 } + l.powi(c_new as i32)
}

/// Processing order (a linear extension) minimizing the modeled DP cost for lattice size l:
/// a shortest path from the empty ideal to the full one in the ideal lattice of p.
fn best_order(p: &Poset, l: f64) -> Vec<usize> {
    let n = p.n;
    let mut ids = ideals(p);
    ids.sort_by_key(|m| m.count_ones());
    let mut best = vec![f64::INFINITY; 1 << n];
    let mut prev = vec![usize::MAX; 1 << n];
    let mut cls = vec![0u32; 1 << n];
    for &m in &ids {
        cls[m as usize] = classes_after(p, m);
    }
    best[0] = 0.0;
    for &m in &ids {
        let m = m as usize;
        for x in 0..n {
            let d = p.down[x] as usize;
            if m >> x & 1 == 1 || d & m != d {
                continue;
            }
            let has_j = (0..n).any(|t| p.cov[t] >> x & 1 == 1);
            let nm = m | 1 << x;
            let c = best[m] + step_cost(l, cls[m], has_j, cls[nm]);
            if c < best[nm] {
                best[nm] = c;
                prev[nm] = x;
            }
        }
    }
    let mut order = Vec::with_capacity(n);
    let mut m = (1usize << n) - 1;
    while m != 0 {
        let x = prev[m];
        order.push(x);
        m ^= 1 << x;
    }
    order.reverse();
    order
}

/// Counter type for the DP: u64 when the total cannot overflow, else u128.
trait Cnt: Copy + PartialEq + std::ops::AddAssign + std::ops::Mul<Output = Self> {
    const ZERO: Self;
    const ONE: Self;
    fn from_usize(x: usize) -> Self;
    fn to_u128(self) -> u128;
}
impl Cnt for u64 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    fn from_usize(x: usize) -> Self {
        x as u64
    }
    fn to_u128(self) -> u128 {
        self as u128
    }
}
impl Cnt for u128 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    fn from_usize(x: usize) -> Self {
        x as u128
    }
    fn to_u128(self) -> u128 {
        self
    }
}

/// The ideal lattice J(P) with the data the Hom counter needs.
struct Lattice {
    masks: Vec<u16>,
    /// idx[mask] = index of the ideal with that mask.
    idx: Vec<u32>,
    /// up[i] = indices of ideals containing ideal i.
    up: Vec<Vec<u32>>,
}

impl Lattice {
    fn new(p: &Poset) -> Self {
        let masks = ideals(p);
        let mut idx = vec![u32::MAX; 1 << p.n];
        for (i, &m) in masks.iter().enumerate() {
            idx[m as usize] = i as u32;
        }
        let up = masks
            .iter()
            .map(|&a| masks.iter().enumerate().filter(|&(_, &b)| a & b == a).map(|(j, _)| j as u32).collect())
            .collect();
        Lattice { masks, idx, up }
    }
    fn len(&self) -> usize {
        self.masks.len()
    }
}

/// Precomputed schedule for the structure side S (elements processed in label order).
///
/// After step j the DP state holds one lattice value per *class*: active elements (assigned,
/// with a later coverer) grouped by their set of future coverers. A future element only needs
/// the join of the values it covers, so a class is represented by the join of its members.
struct Step {
    /// old classes whose join feeds element j's lower bound
    feeds: Vec<usize>,
    /// for each new class: the old classes merged into it, and whether j itself belongs to it
    new_classes: Vec<(Vec<usize>, bool)>,
}

struct Plan {
    steps: Vec<Step>,
    max_classes: usize,
    /// (classes before, j keeps a class, classes after) per step, for the cost model
    shape: Vec<(u32, bool, u32)>,
}

impl Plan {
    fn new(s: &Poset) -> Self {
        let k = s.n;
        let cov = &s.cov[..k];
        // future coverers of i after step j
        let fut = |i: usize, j: usize| -> u16 { (j + 1..k).filter(|&t| cov[t] >> i & 1 == 1).fold(0, |m, t| m | 1 << t) };
        let mut classes: Vec<u16> = Vec::new(); // future-coverer set per class (after previous step)
        let mut steps = Vec::with_capacity(k);
        let mut max_classes = 0;
        let mut shape = Vec::with_capacity(k);
        for j in 0..k {
            let c_old = classes.len() as u32;
            let feeds: Vec<usize> = (0..classes.len()).filter(|&c| classes[c] >> j & 1 == 1).collect();
            let mut next: Vec<(u16, Vec<usize>, bool)> = Vec::new();
            let put = |set: u16, old: Option<usize>, next: &mut Vec<(u16, Vec<usize>, bool)>| {
                if set == 0 {
                    return;
                }
                let pos = match next.iter().position(|e| e.0 == set) {
                    Some(p) => p,
                    None => {
                        next.push((set, Vec::new(), false));
                        next.len() - 1
                    }
                };
                match old {
                    Some(c) => next[pos].1.push(c),
                    None => next[pos].2 = true,
                }
            };
            for (c, &set) in classes.iter().enumerate() {
                put(set & !(1 << j), Some(c), &mut next);
            }
            put(fut(j, j), None, &mut next);
            classes = next.iter().map(|e| e.0).collect();
            max_classes = max_classes.max(classes.len());
            shape.push((c_old, next.iter().any(|e| e.2), classes.len() as u32));
            steps.push(Step { feeds, new_classes: next.into_iter().map(|e| (e.1, e.2)).collect() });
        }
        assert!(classes.is_empty());
        Plan { steps, max_classes, shape }
    }
    /// Plan for p, processed in the order that is cheapest over lattices of size about l.
    fn optimized(p: &Poset, l: f64) -> Self {
        Plan::new(&relabel(p, &best_order(p, l)))
    }
    fn cost(&self, l: usize) -> f64 {
        self.shape.iter().map(|&(a, h, b)| step_cost(l as f64, a, h, b)).sum()
    }
}

/// Hom(S, L): number of order-preserving maps S -> L, by a DP over S in label order whose state
/// is one lattice value (ideal index) per class, stored dense in mixed radix |L|^classes.
fn hom<T: Cnt>(plan: &Plan, lat: &Lattice, buf: &mut [Vec<T>; 2]) -> u128 {
    let l = lat.len();
    let mut nold = 0usize;
    buf[0].clear();
    buf[0].push(T::ONE);
    let mut vals = [0usize; MAXN];
    for st in &plan.steps {
        let nnew = st.new_classes.len();
        let (a, b) = buf.split_at_mut(1);
        let (old, new) = (&a[0], &mut b[0]);
        new.clear();
        new.resize(l.pow(nnew as u32), T::ZERO);
        let j_class = st.new_classes.iter().position(|c| c.1);
        for (si, &c) in old.iter().enumerate() {
            if c == T::ZERO {
                continue;
            }
            let mut r = si;
            for v in vals.iter_mut().take(nold) {
                *v = r % l;
                r /= l;
            }
            let mut lo = 0u16;
            for &f in &st.feeds {
                lo |= lat.masks[vals[f]];
            }
            let ups = &lat.up[lat.idx[lo as usize] as usize];
            if j_class.is_none() {
                let mut base = 0usize;
                let mut mul = 1usize;
                for (members, _) in &st.new_classes {
                    let m = members.iter().fold(0u16, |m, &o| m | lat.masks[vals[o]]);
                    base += lat.idx[m as usize] as usize * mul;
                    mul *= l;
                }
                new[base] += c * T::from_usize(ups.len());
            } else {
                // j's class value is join(j's value, merged old members): compute the rest once
                let mut base = 0usize;
                let mut mul = 1usize;
                let mut jmul = 0usize;
                let mut jrest = 0u16;
                for (members, has_j) in &st.new_classes {
                    let m = members.iter().fold(0u16, |m, &o| m | lat.masks[vals[o]]);
                    if *has_j {
                        jmul = mul;
                        jrest = m;
                    } else {
                        base += lat.idx[m as usize] as usize * mul;
                    }
                    mul *= l;
                }
                if jrest == 0 {
                    for &u in ups {
                        new[base + u as usize * jmul] += c;
                    }
                } else {
                    for &u in ups {
                        let m = lat.masks[u as usize] | jrest;
                        new[base + lat.idx[m as usize] as usize * jmul] += c;
                    }
                }
            }
        }
        buf.swap(0, 1);
        nold = nnew;
    }
    buf[0][0].to_u128()
}

/// Print the matrix M[P][Q] = Hom(Q, J(P)) over unlabeled posets P on n and Q on k points, with
/// the weights, as TSV: first line the Q weights, then one line per P (weight, then the row).
/// For exploring structure (e.g. the rank of M).
fn cmd_matrix(n: usize, k: usize) {
    let ps = load(n);
    let qs = load(k);
    let plans: Vec<Plan> = qs.iter().map(Plan::new).collect();
    let mut buf: [Vec<u128>; 2] = [Vec::new(), Vec::new()];
    let qw: Vec<String> = qs.iter().map(|q| weight(q).to_string()).collect();
    println!("w\t{}", qw.join("\t"));
    for p in &ps {
        let lat = Lattice::new(p);
        let row: Vec<String> = plans.iter().map(|pl| hom(pl, &lat, &mut buf).to_string()).collect();
        println!("{}\t{}", weight(p), row.join("\t"));
    }
}

fn cmd_check(nmax: usize) {
    for n in 1..=nmax {
        let t = Instant::now();
        let ps = load(n);
        let fact: u128 = (1..=n as u128).product();
        let (mut sw, mut sl) = (0u128, 0u128);
        for p in &ps {
            sw += weight(p);
            sl += fact / automorphisms(p);
        }
        let ok_w = n >= A006455.len() || sw == A006455[n];
        let ok_l = n >= A001035.len() || sl == A001035[n];
        println!(
            "n={n}: {} posets, sum w = {sw} [{}], sum n!/|Aut| = {sl} [{}]  ({:.2}s)",
            ps.len(),
            if ok_w { "ok" } else { "MISMATCH" },
            if ok_l { "ok" } else { "MISMATCH" },
            t.elapsed().as_secs_f64()
        );
        assert!(ok_w && ok_l);
    }
}

fn cmd_split(n: usize, k: usize, threads: usize, flags: &str) {
    let t0 = Instant::now();
    let (opt_order, opt_dir, opt_sym) = (flags.contains('o'), flags.contains('d'), flags.contains('s'));
    let ps = load(n);
    let ss = load(k);
    let pw: Vec<u128> = ps.iter().map(weight).collect();
    // lattice sizes for the order heuristic: the mean ideal count of the other side
    let mean_ideals = |v: &[Poset]| v.iter().map(|p| ideals(p).len() as f64).sum::<f64>() / v.len() as f64;
    let (lp, ls) = (mean_ideals(&ps), mean_ideals(&ss));
    let plan = |p: &Poset, l: f64| if opt_order { Plan::optimized(p, l) } else { Plan::new(p) };
    // S side: plan for S over J(P), and the lattice J(S^op) for the reverse direction
    let sdata: Vec<(u128, Plan, Lattice)> =
        ss.iter().map(|s| (weight(s), plan(s, lp), Lattice::new(&dual(s)))).collect();
    let maxf = sdata.iter().map(|d| d.1.max_classes).max().unwrap();
    let wide = n * k >= 64; // Hom <= 2^(nk): u128 counters needed
    let sym = opt_sym && n == k; // the summand is symmetric in (P, S)
    eprintln!(
        "split {n}+{k} [{flags}]: {} x {} posets, max classes {maxf}, {} counters{}, setup {:.2}s",
        ps.len(),
        ss.len(),
        if wide { "u128" } else { "u64" },
        if sym { ", unordered pairs" } else { "" },
        t0.elapsed().as_secs_f64()
    );
    let next = AtomicUsize::new(0);
    let flipped = AtomicUsize::new(0);
    let total: u128 = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                sc.spawn(|| {
                    let mut b64: [Vec<u64>; 2] = [Vec::new(), Vec::new()];
                    let mut b128: [Vec<u128>; 2] = [Vec::new(), Vec::new()];
                    let mut acc = 0u128;
                    let mut nflip = 0usize;
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= ps.len() {
                            break;
                        }
                        let lat = Lattice::new(&ps[i]);
                        // P^op over J(S^op): the same count, i(S^op x P) = i(P x S^op)
                        let pplan = if opt_dir { Some(plan(&dual(&ps[i]), ls)) } else { None };
                        let mut inner = 0u128;
                        let range = if sym { i..ss.len() } else { 0..ss.len() };
                        for jdx in range {
                            let (w, splan, slat) = &sdata[jdx];
                            let (pl, l) = match &pplan {
                                Some(pp) if pp.cost(slat.len()) < splan.cost(lat.len()) => {
                                    nflip += 1;
                                    (pp, slat)
                                }
                                _ => (splan, &lat),
                            };
                            let h = if wide { hom(pl, l, &mut b128) } else { hom(pl, l, &mut b64) };
                            let mult = if sym && jdx != i { 2 } else { 1 };
                            inner += mult * w * h;
                        }
                        acc += pw[i] * inner;
                    }
                    flipped.fetch_add(nflip, Ordering::Relaxed);
                    acc
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).sum()
    });
    let nn = n + k;
    let status = if nn < A006455.len() {
        if total == A006455[nn] { "matches OEIS" } else { "MISMATCH with OEIS" }
    } else {
        "new"
    };
    eprintln!("  {} pairs ran in the reverse direction", flipped.load(Ordering::Relaxed));
    println!("a({nn}) = {total}  [split {n}+{k} {flags}, {status}, {:.2}s]", t0.elapsed().as_secs_f64());
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
    let args: Vec<String> = std::env::args().collect();
    let num = |i: usize, d: usize| args.get(i).map_or(d, |s| s.parse().expect("number"));
    match args.get(1).map(String::as_str) {
        Some("check") => cmd_check(num(2, 9)),
        Some("matrix") => cmd_matrix(num(2, 0), num(3, 0)),
        Some("split") => cmd_split(num(2, 0), num(3, 0), num(4, 10), args.get(5).map_or("ods", |s| s.as_str())),
        _ => eprintln!("usage: a006455 check [nmax] | split n k [threads] [flags]"),
    }
}
