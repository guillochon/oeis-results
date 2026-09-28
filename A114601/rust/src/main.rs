//! A114601: n x n symmetric positive definite matrices with 2 on the diagonal and entries in
//! {-1, 0, 1} elsewhere. Two independent routes to the E_8 term of the root-lattice formula
//! (see ../README.md):
//!
//!   a114601 gram <n> [--ckpt FILE]
//!                            Route 1: depth-first enumeration of all PD matrices of size n.
//!                            Prints a(n) and the connected matrices by determinant (E_8: det 1,
//!                            D_8: det 4, A_8: det 9 at n = 8).
//!   a114601 roots <E6|E7|E8> [--fix]
//!                            Route 2: n-subsets of the positive roots of E_n (actual vectors),
//!                            counted by Gram determinant; bases have det = det(E_n). --fix keeps
//!                            only subsets containing root 0 (W is transitive on roots), which is
//!                            multiplied back by (#positive roots)/n. Prints the E-term
//!                            #bases * 2^n * n! / |Aut(E_n)| and a Cauchy-Binet check.
//!
//! Exact arithmetic throughout: with A = adj(M_k) and d = det(M_k), appending a row r with diagonal
//! 2 gives det(M_{k+1}) = 2d - r^T A r, and
//!     adj(M_{k+1}) = [[(det' A + u u^T) / d, -u], [-u^T, d]],   u = A r.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const NMAX: usize = 9;
const DMAX: usize = 1 << NMAX; // determinants of these Gram matrices are at most 2^n

#[derive(Clone)]
struct Node {
    k: usize,
    m: [[i8; NMAX]; NMAX],
    det: i64,
    adj: [[i64; NMAX]; NMAX],
}

impl Node {
    fn root() -> Node {
        let mut n = Node { k: 1, m: [[0; NMAX]; NMAX], det: 2, adj: [[0; NMAX]; NMAX] };
        n.m[0][0] = 2;
        n.adj[0][0] = 1;
        n
    }

    /// r^T A r
    #[inline]
    fn quad(&self, r: &[i8]) -> i64 {
        let k = self.k;
        let mut q = 0i64;
        for i in 0..k {
            if r[i] == 0 {
                continue;
            }
            let mut s = 0i64;
            for j in 0..k {
                s += self.adj[i][j] * r[j] as i64;
            }
            q += s * r[i] as i64;
        }
        q
    }

    /// The child obtained by appending row r (with new determinant nd > 0).
    fn child(&self, r: &[i8], nd: i64) -> Node {
        let k = self.k;
        let mut c = self.clone();
        let mut u = [0i64; NMAX];
        for i in 0..k {
            u[i] = (0..k).map(|j| self.adj[i][j] * r[j] as i64).sum();
        }
        for i in 0..k {
            for j in 0..k {
                let v = nd * self.adj[i][j] + u[i] * u[j];
                debug_assert_eq!(v % self.det, 0);
                c.adj[i][j] = v / self.det;
            }
            c.adj[i][k] = -u[i];
            c.adj[k][i] = -u[i];
            c.m[i][k] = r[i];
            c.m[k][i] = r[i];
        }
        c.adj[k][k] = self.det;
        c.m[k][k] = 2;
        c.det = nd;
        c.k = k + 1;
        c
    }

    /// Bitmasks of the connected components of the graph of nonzero off-diagonal entries.
    fn components(&self) -> Vec<u32> {
        let k = self.k;
        let mut seen = 0u32;
        let mut out = vec![];
        for s in 0..k {
            if seen >> s & 1 == 1 {
                continue;
            }
            let mut comp = 1u32 << s;
            let mut stack = vec![s];
            while let Some(i) = stack.pop() {
                for j in 0..k {
                    if self.m[i][j] != 0 && comp >> j & 1 == 0 {
                        comp |= 1 << j;
                        stack.push(j);
                    }
                }
            }
            seen |= comp;
            out.push(comp);
        }
        out
    }
}

fn candidates(k: usize) -> Vec<[i8; NMAX]> {
    let mut out = vec![];
    let total = 3usize.pow(k as u32);
    for mut x in 0..total {
        let mut r = [0i8; NMAX];
        for slot in r.iter_mut().take(k) {
            *slot = (x % 3) as i8 - 1;
            x /= 3;
        }
        out.push(r);
    }
    out
}

// ---------------------------------------------------------------------------------------------
// Route 1: all PD Gram matrices of size n
// ---------------------------------------------------------------------------------------------

struct Counts {
    total: u64,
    conn: Vec<u64>, // connected matrices by determinant
}

/// Last level: enumerate the final row recursively, accumulating r^T A r and the support mask.
fn leaves(node: &Node, comps: &[u32], cnt: &mut Counts) {
    let k = node.k;
    fn rec(node: &Node, comps: &[u32], i: usize, p: &mut [i64; NMAX], q: i64, mask: u32, cnt: &mut Counts) {
        let k = node.k;
        if i == k {
            let nd = 2 * node.det - q;
            if nd > 0 {
                cnt.total += 1;
                if comps.iter().all(|&c| c & mask != 0) {
                    cnt.conn[nd as usize] += 1;
                }
            }
            return;
        }
        // r_i = 0
        rec(node, comps, i + 1, p, q, mask, cnt);
        for s in [-1i64, 1] {
            // q += r_i^2 A_ii + 2 r_i sum_{j<i} r_j A_ij  (p[i] = sum_{j<i} r_j A_ji)
            let nq = q + node.adj[i][i] + 2 * s * p[i];
            for j in i + 1..k {
                p[j] += s * node.adj[i][j];
            }
            rec(node, comps, i + 1, p, nq, mask | 1 << i, cnt);
            for j in i + 1..k {
                p[j] -= s * node.adj[i][j];
            }
        }
    }
    let mut p = [0i64; NMAX];
    let _ = k;
    rec(node, comps, 0, &mut p, 0, 0, cnt);
}

fn gram_subtree(node: &Node, n: usize, cands: &[Vec<[i8; NMAX]>], cnt: &mut Counts) {
    if node.k == n {
        cnt.total += 1;
        if node.components().len() == 1 {
            cnt.conn[node.det as usize] += 1;
        }
        return;
    }
    if node.k == n - 1 {
        let comps = node.components();
        leaves(node, &comps, cnt);
        return;
    }
    for r in &cands[node.k] {
        let nd = 2 * node.det - node.quad(r);
        if nd > 0 {
            gram_subtree(&node.child(r, nd), n, cands, cnt);
        }
    }
}

fn gram(n: usize) {
    assert!((1..=NMAX).contains(&n));
    let t0 = Instant::now();
    let cands: Vec<Vec<[i8; NMAX]>> = (0..NMAX).map(candidates).collect();
    // split: all nodes at level L, then subtrees in parallel
    let level = if n >= 6 { n - 3 } else { 1 };
    let mut frontier = vec![Node::root()];
    while frontier[0].k < level {
        let mut next = vec![];
        for nd in &frontier {
            for r in &cands[nd.k] {
                let d = 2 * nd.det - nd.quad(r);
                if d > 0 {
                    next.push(nd.child(r, d));
                }
            }
        }
        frontier = next;
    }
    // Subtrees are processed in chunks of CHUNK frontier nodes. With --ckpt FILE, every finished
    // chunk appends "chunk total det:count ..." to FILE, and a rerun skips those chunks.
    const CHUNK: usize = 1000;
    let nchunks = frontier.len().div_ceil(CHUNK);
    let ckpt = CKPT.lock().unwrap().clone();
    let total = Mutex::new(Counts { total: 0, conn: vec![0; DMAX + 1] });
    let mut done = vec![false; nchunks];
    if let Some(path) = &ckpt {
        if let Ok(text) = std::fs::read_to_string(path) {
            let mut t = total.lock().unwrap();
            for line in text.lines() {
                let f: Vec<&str> = line.split_whitespace().collect();
                if f.len() < 2 || !line.ends_with(';') {
                    continue; // incomplete line
                }
                let ch: usize = f[0].parse().unwrap();
                if ch >= nchunks || done[ch] {
                    continue;
                }
                done[ch] = true;
                t.total += f[1].parse::<u64>().unwrap();
                for kv in &f[2..] {
                    let kv = kv.trim_end_matches(';');
                    if let Some((d, c)) = kv.split_once(':') {
                        t.conn[d.parse::<usize>().unwrap()] += c.parse::<u64>().unwrap();
                    }
                }
            }
            eprintln!("  resuming: {} of {nchunks} chunks already in {path}", done.iter().filter(|&&x| x).count());
        }
    }
    let todo: Vec<usize> = (0..nchunks).filter(|&c| !done[c]).collect();
    let writer = ckpt.as_ref().map(|p| Mutex::new(std::fs::OpenOptions::new().create(true).append(true).open(p).unwrap()));
    let idx = AtomicUsize::new(0);
    let finished = AtomicUsize::new(nchunks - todo.len());
    std::thread::scope(|sc| {
        for _ in 0..num_threads() {
            sc.spawn(|| loop {
                let ti = idx.fetch_add(1, Ordering::Relaxed);
                if ti >= todo.len() {
                    break;
                }
                let ch = todo[ti];
                let mut c = Counts { total: 0, conn: vec![0; DMAX + 1] };
                for node in &frontier[ch * CHUNK..((ch + 1) * CHUNK).min(frontier.len())] {
                    gram_subtree(node, n, &cands, &mut c);
                }
                if let Some(w) = &writer {
                    use std::io::Write;
                    let mut line = format!("{ch} {}", c.total);
                    for (d, &x) in c.conn.iter().enumerate() {
                        if x > 0 {
                            line += &format!(" {d}:{x}");
                        }
                    }
                    line += ";\n";
                    let mut f = w.lock().unwrap();
                    f.write_all(line.as_bytes()).unwrap();
                    f.flush().unwrap();
                }
                let fin = finished.fetch_add(1, Ordering::Relaxed) + 1;
                if n >= 8 {
                    let el = t0.elapsed().as_secs_f64();
                    let run = fin - (nchunks - todo.len());
                    eprintln!(
                        "  [{el:7.0}s] chunks {fin}/{nchunks}  (ETA {:.1} h)",
                        el / run as f64 * (nchunks - fin) as f64 / 3600.0
                    );
                }
                let mut t = total.lock().unwrap();
                t.total += c.total;
                for (a, b) in t.conn.iter_mut().zip(&c.conn) {
                    *a += b;
                }
            });
        }
    });
    let t = total.into_inner().unwrap();
    let conn: Vec<(usize, u64)> = t.conn.iter().enumerate().filter(|(_, &c)| c > 0).map(|(d, &c)| (d, c)).collect();
    println!(
        "gram n={n}: a({n}) = {}   connected by det: {:?}   connected total {}   ({:.1}s)",
        t.total,
        conn,
        conn.iter().map(|x| x.1).sum::<u64>(),
        t0.elapsed().as_secs_f64()
    );
}

// ---------------------------------------------------------------------------------------------
// Route 2: bases of E_6, E_7, E_8 made of roots
// ---------------------------------------------------------------------------------------------

/// The 240 roots of E_8, in doubled coordinates (so entries are integers, norm 8).
fn e8_roots() -> Vec<[i32; 8]> {
    let mut out = vec![];
    for i in 0..8 {
        for j in i + 1..8 {
            for si in [2, -2] {
                for sj in [2, -2] {
                    let mut v = [0; 8];
                    v[i] = si;
                    v[j] = sj;
                    out.push(v);
                }
            }
        }
    }
    for m in 0..256u32 {
        if m.count_ones() % 2 == 0 {
            let mut v = [1; 8];
            for (i, x) in v.iter_mut().enumerate() {
                if m >> i & 1 == 1 {
                    *x = -1;
                }
            }
            out.push(v);
        }
    }
    assert_eq!(out.len(), 240);
    out
}

fn ip(a: &[i32; 8], b: &[i32; 8]) -> i32 {
    a.iter().zip(b).map(|(x, y)| x * y).sum::<i32>() / 4
}

fn roots(kind: &str, fix: bool) {
    let t0 = Instant::now();
    let all = e8_roots();
    // E_7 = roots orthogonal to one root; E_6 = roots orthogonal to an A_2 (two roots with ip -1)
    let alpha = all[0];
    let beta = *all.iter().find(|r| ip(r, &alpha) == -1).unwrap();
    let (n, sel, det_full, aut, h): (usize, Vec<[i32; 8]>, i64, u128, i64) = match kind {
        "E8" => (8, all.clone(), 1, 696_729_600, 30),
        "E7" => (7, all.iter().copied().filter(|r| ip(r, &alpha) == 0).collect(), 2, 2_903_040, 18),
        "E6" => (6, all.iter().copied().filter(|r| ip(r, &alpha) == 0 && ip(r, &beta) == 0).collect(), 3, 103_680, 12),
        _ => panic!("kind must be E6, E7 or E8"),
    };
    // one representative per +- pair
    let pos: Vec<[i32; 8]> = sel.iter().copied().filter(|v| v.iter().find(|&&x| x != 0).unwrap() > &0).collect();
    let np = pos.len();
    let g: Vec<Vec<i8>> = pos.iter().map(|a| pos.iter().map(|b| ip(a, b) as i8).collect()).collect();
    println!("{kind}: {} roots, {np} positive; subsets of size {n}{}", sel.len(), if fix { " containing root 0" } else { "" });

    // DFS over increasing index sequences; tasks = first two indices
    let mut tasks = vec![];
    let first: Vec<usize> = if fix { vec![0] } else { (0..np).collect() };
    for &a in &first {
        for b in a + 1..np {
            tasks.push((a, b));
        }
    }
    let idx = AtomicUsize::new(0);
    let hist = Mutex::new(vec![0u64; 4096]);
    std::thread::scope(|sc| {
        for _ in 0..num_threads() {
            sc.spawn(|| {
                let mut h = vec![0u64; 4096];
                loop {
                    let i = idx.fetch_add(1, Ordering::Relaxed);
                    if i >= tasks.len() {
                        break;
                    }
                    let (a, b) = tasks[i];
                    let mut node = Node::root();
                    let r = [g[a][b]; 1];
                    let mut rr = [0i8; NMAX];
                    rr[0] = r[0];
                    let nd = 2 * node.det - node.quad(&rr);
                    if nd <= 0 {
                        continue;
                    }
                    node = node.child(&rr, nd);
                    let mut chosen = vec![a, b];
                    subsets(&g, &mut chosen, &node, n, np, &mut h);
                }
                let mut t = hist.lock().unwrap();
                for (x, y) in t.iter_mut().zip(&h) {
                    *x += y;
                }
            });
        }
    });
    let hist = hist.into_inner().unwrap();
    let mut by_det: Vec<(usize, u64)> = hist.iter().enumerate().filter(|(_, &c)| c > 0).map(|(d, &c)| (d, c)).collect();
    by_det.sort();
    let sum_det: u128 = by_det.iter().map(|&(d, c)| d as u128 * c as u128).sum();
    let bases = hist[det_full as usize] as u128;
    // scale back from "subsets containing root 0" to all subsets
    let (bases_all, sum_all) = if fix {
        assert_eq!((bases * np as u128) % n as u128, 0);
        (bases * np as u128 / n as u128, sum_det * np as u128 / n as u128)
    } else {
        (bases, sum_det)
    };
    let cb = (h as u128).pow(n as u32);
    let ordered = bases_all * (1u128 << n) * (1..=n as u128).product::<u128>();
    println!("  independent subsets by Gram det: {by_det:?}");
    println!("  Cauchy-Binet: sum of dets = {sum_all}, expected h^n = {h}^{n} = {cb}  {}", if sum_all == cb { "OK" } else { "MISMATCH" });
    println!(
        "  bases (det {det_full}): {bases_all} unordered sets of positive roots -> {ordered} ordered root bases; / |Aut| {aut} = {} (remainder {})",
        ordered / aut,
        ordered % aut
    );
    println!("  => c_{kind} = {}   ({:.1}s)", ordered / aut, t0.elapsed().as_secs_f64());
}

fn subsets(g: &[Vec<i8>], chosen: &mut Vec<usize>, node: &Node, n: usize, np: usize, h: &mut [u64]) {
    if node.k == n {
        h[node.det as usize] += 1;
        return;
    }
    let last = *chosen.last().unwrap();
    let mut r = [0i8; NMAX];
    for c in last + 1..np {
        for (t, &x) in chosen.iter().enumerate() {
            r[t] = g[x][c];
        }
        let nd = 2 * node.det - node.quad(&r);
        if nd <= 0 {
            continue; // dependent (Gram of real vectors: det >= 0, 0 iff dependent)
        }
        if node.k + 1 == n {
            h[nd as usize] += 1;
        } else {
            chosen.push(c);
            subsets(g, chosen, &node.child(&r, nd), n, np, h);
            chosen.pop();
        }
    }
}

// ---------------------------------------------------------------------------------------------

static THREADS: AtomicUsize = AtomicUsize::new(0);
static CKPT: Mutex<Option<String>> = Mutex::new(None);

fn num_threads() -> usize {
    match THREADS.load(Ordering::Relaxed) {
        0 => std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).saturating_sub(2).max(1),
        t => t,
    }
}

/// Windows 11 keeps background console jobs on the E-cores (EcoQoS); opt out.
#[cfg(windows)]
fn disable_power_throttling() {
    #[repr(C)]
    struct State {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
    }
    let st = State { version: 1, control_mask: 1, state_mask: 0 };
    unsafe {
        SetProcessInformation(GetCurrentProcess(), 4, &st as *const _ as *const std::ffi::c_void, std::mem::size_of::<State>() as u32);
    }
}

#[cfg(not(windows))]
fn disable_power_throttling() {}

fn main() {
    disable_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Some(i) = args.iter().position(|a| a == "--threads") {
        THREADS.store(args[i + 1].parse().unwrap(), Ordering::Relaxed);
    }
    if let Some(i) = args.iter().position(|a| a == "--ckpt") {
        *CKPT.lock().unwrap() = Some(args[i + 1].clone());
    }
    match args.first().map(String::as_str) {
        Some("gram") => gram(args[1].parse().unwrap()),
        Some("roots") => roots(&args[1], args.iter().any(|a| a == "--fix")),
        _ => eprintln!("usage: a114601 gram <n> | roots <E6|E7|E8> [--fix]"),
    }
}
