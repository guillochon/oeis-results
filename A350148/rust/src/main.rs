//! Frontier transfer matrix for A350148 (Hilbert-style square-filling motifs, McKenna).
//!
//! Input: an n x n grid of frozen perimeter tiles ('B','R','T','L' = which edge the tile owns, '.' = free),
//! written by `make_bounds.py` for one Fibbinary zipping mode. Counts the ways to give every free tile one of its
//! four edges so that the n^2 owned edges are pairwise distinct and form a path from a = (0,0) to b = (n,0) that
//! visits every interior lattice point (toroidal Hamiltonicity, McKenna Lemma 3.17). Boundary lattice points are
//! touched only by perimeter tiles, so the frozen boundary fixes which of them the path uses.
//!
//! The a-b path plus a virtual edge b-a is a cycle; the DP is the usual Hamiltonian-cycle frontier DP with component
//! labels, plus one "claimed" bit per frontier edge (an edge owned by the tile on one side cannot be owned by the
//! other) and a dead-end check when a lattice point leaves the frontier.
//!
//! Modes:  count FILE            -> number of motifs
//!         enum FILE             -> every motif, one line of n^2 edge letters (row-major from the bottom)
//!         sample FILE K SEED    -> K uniformly random motifs

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
struct Mix(u64);
impl Hasher for Mix {
    fn finish(&self) -> u64 {
        let mut z = self.0.wrapping_add(0x9E3779B97F4A7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ b as u64;
        }
    }
    fn write_u128(&mut self, v: u128) {
        self.0 = (v as u64) ^ ((v >> 64) as u64).wrapping_mul(0x9E3779B97F4A7C15);
    }
}
type Map = HashMap<u128, u128, BuildHasherDefault<Mix>>;

const MAXP: usize = 24; // n + 1 <= MAXP
const DEG2: u8 = 15; // vertex code: 0 = unused so far, 1..14 = one end of strand with that label, 15 = degree 2

/// Frontier: p[0..=n] as described in the module doc, k = the vertex above-left of the next cell,
/// h[i] = frontier horizontal edge i already owned, v = vertical edge left of the next cell already owned.
#[derive(Clone)]
struct St {
    p: [u8; MAXP],
    k: u8,
    h: u32,
    v: bool,
    n: u8,
}

fn encode(s: &St) -> u128 {
    let mut x: u128 = 0;
    for &c in s.p[..=s.n as usize].iter().chain(std::iter::once(&s.k)) {
        x = (x << 4) | c as u128;
    }
    x = (x << 32) | s.h as u128;
    (x << 1) | s.v as u128
}

fn decode(mut x: u128, n: usize) -> St {
    let v = x & 1 == 1;
    x >>= 1;
    let h = (x & 0xffff_ffff) as u32;
    x >>= 32;
    let k = (x & 15) as u8;
    x >>= 4;
    let mut p = [0u8; MAXP];
    for i in (0..=n).rev() {
        p[i] = (x & 15) as u8;
        x >>= 4;
    }
    St { p, k, h, v, n: n as u8 }
}

/// Relabel strand ends in order of first appearance.
fn canon(s: &mut St) {
    let mut map = [0u8; 16];
    let mut next = 1u8;
    let n = s.n as usize;
    for c in s.p[..=n].iter_mut().chain(std::iter::once(&mut s.k)) {
        if *c != 0 && *c != DEG2 {
            if map[*c as usize] == 0 {
                map[*c as usize] = next;
                next += 1;
            }
            *c = map[*c as usize];
        }
    }
}

/// Vertex slots used while processing a cell: 0 = BL (p[x]), 1 = BR (p[x+1]), 2 = TL (k), 3 = TR (new).
struct Work {
    c: [u8; 4],
    rest: [u8; MAXP], // a copy of p, to relabel the far ends of strands
}

fn label_of(w: &Work, i: usize) -> u8 {
    w.c[i]
}

/// Add edge between slots i and j. Returns false if it makes a degree 3 or closes a cycle (unless allowed).
fn add_edge(w: &mut Work, x: usize, i: usize, j: usize, may_close: bool) -> bool {
    let (a, b) = (label_of(w, i), label_of(w, j));
    if a == DEG2 || b == DEG2 {
        return false;
    }
    match (a == 0, b == 0) {
        (true, true) => {
            let l = 14; // fresh label, canonicalised later
            w.c[i] = l;
            w.c[j] = l;
        }
        (true, false) => {
            w.c[i] = b;
            w.c[j] = DEG2;
        }
        (false, true) => {
            w.c[j] = a;
            w.c[i] = DEG2;
        }
        (false, false) => {
            if a == b {
                if !may_close {
                    return false;
                }
                w.c[i] = DEG2;
                w.c[j] = DEG2;
            } else {
                w.c[i] = DEG2;
                w.c[j] = DEG2;
                // the other end of strand b becomes strand a
                for q in 0..4 {
                    if w.c[q] == b {
                        w.c[q] = a;
                    }
                }
                for (idx, c) in w.rest.iter_mut().enumerate() {
                    if idx != x && idx != x + 1 && *c == b {
                        *c = a;
                    }
                }
            }
        }
    }
    true
}

struct Problem {
    n: usize,
    frozen: Vec<Option<u8>>, // per cell (row-major from bottom): forced edge 0=B 1=R 2=T 3=L
}

fn interior(n: usize, i: usize, j: usize) -> bool {
    i > 0 && i < n && j > 0 && j < n
}

/// Leaving lattice point with code c: interior points must have degree 2; boundary points 0 or 2.
fn final_ok(n: usize, i: usize, j: usize, c: u8) -> bool {
    if interior(n, i, j) {
        c == DEG2
    } else {
        c == 0 || c == DEG2
    }
}

/// All successors of state s when processing cell idx. Each successor carries the edge chosen (0..3).
fn step(pr: &Problem, idx: usize, s: &St, out: &mut Vec<(u8, St)>) {
    let n = pr.n;
    let (x, y) = (idx % n, idx / n);
    let last = idx == n * n - 1;
    let choices: &[u8] = match pr.frozen[idx] {
        Some(e) => match e {
            0 => &[0],
            1 => &[1],
            2 => &[2],
            _ => &[3],
        },
        None => &[0, 1, 2, 3],
    };
    for &e in choices {
        let hx = s.h >> x & 1 == 1;
        if e == 0 && hx {
            continue;
        }
        if e == 3 && s.v {
            continue;
        }
        let mut w = Work { c: [s.p[x], s.p[x + 1], s.k, 0], rest: s.p };
        let ok = match e {
            0 => add_edge(&mut w, x, 0, 1, last),
            1 => add_edge(&mut w, x, 1, 3, last),
            2 => add_edge(&mut w, x, 2, 3, last),
            _ => add_edge(&mut w, x, 0, 2, last),
        };
        if !ok {
            continue;
        }
        // BL = (x, y) leaves
        if !final_ok(n, x, y, w.c[0]) {
            continue;
        }
        let mut p = w.rest;
        p[x] = w.c[2];
        p[x + 1] = w.c[1];
        let mut h = s.h & !(1 << x);
        if e == 2 {
            h |= 1 << x;
        }
        let mut ns = St { p, k: w.c[3], h, v: e == 1, n: n as u8 };
        if x == n - 1 {
            // (n, y) leaves; the new row starts with a fresh top-left point
            if !final_ok(n, n, y, ns.p[n]) {
                continue;
            }
            ns.p[n] = ns.k;
            ns.k = 0;
            ns.v = false;
            if y == n - 1 {
                // whole top line leaves
                if !(0..=n).all(|i| final_ok(n, i, n, ns.p[i])) {
                    continue;
                }
            }
        }
        canon(&mut ns);
        out.push((e, ns));
    }
}

fn initial(n: usize) -> St {
    let mut p = [0u8; MAXP];
    p[0] = 1; // a and b: one end each of the strand closed by the virtual edge b-a
    p[n] = 1;
    St { p, k: 0, h: 0, v: false, n: n as u8 }
}

fn read_problem(path: &str) -> Problem {
    let txt = std::fs::read_to_string(path).expect("read bounds file");
    let rows: Vec<&str> = txt.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')).collect();
    let n = rows.len();
    let mut frozen = vec![None; n * n];
    // file lists the top row first
    for (r, line) in rows.iter().enumerate() {
        let y = n - 1 - r;
        for (x, ch) in line.trim().chars().enumerate() {
            frozen[y * n + x] = match ch {
                'B' => Some(0),
                'R' => Some(1),
                'T' => Some(2),
                'L' => Some(3),
                _ => None,
            };
        }
    }
    Problem { n, frozen }
}

/// Forward pass. If keep, returns every layer's map (layer i = states before cell i).
fn forward(pr: &Problem, keep: bool) -> (u128, Vec<Map>) {
    let n = pr.n;
    let mut cur: Map = Map::default();
    cur.insert(encode(&initial(n)), 1);
    let mut layers = Vec::new();
    let mut buf = Vec::new();
    let mut peak = 0usize;
    for idx in 0..n * n {
        let mut nxt: Map = Map::default();
        for (&code, &cnt) in cur.iter() {
            let s = decode(code, n);
            buf.clear();
            step(pr, idx, &s, &mut buf);
            for (_, ns) in buf.drain(..) {
                *nxt.entry(encode(&ns)).or_insert(0) += cnt;
            }
        }
        peak = peak.max(nxt.len());
        if keep {
            layers.push(std::mem::replace(&mut cur, nxt));
        } else {
            cur = nxt;
        }
    }
    eprintln!("peak states {}", peak);
    let total: u128 = cur.values().sum();
    if keep {
        layers.push(cur);
    }
    (total, layers)
}

/// ways[i][s] = number of completions from state s before cell i.
fn backward(pr: &Problem, layers: &[Map]) -> Vec<Map> {
    let n = pr.n;
    let mut ways: Vec<Map> = vec![Map::default(); n * n + 1];
    for (&code, _) in layers[n * n].iter() {
        ways[n * n].insert(code, 1);
    }
    let mut buf = Vec::new();
    for idx in (0..n * n).rev() {
        let mut w = Map::default();
        for (&code, _) in layers[idx].iter() {
            let s = decode(code, n);
            buf.clear();
            step(pr, idx, &s, &mut buf);
            let mut tot = 0u128;
            for (_, ns) in buf.drain(..) {
                tot += ways[idx + 1].get(&encode(&ns)).copied().unwrap_or(0);
            }
            if tot > 0 {
                w.insert(code, tot);
            }
        }
        ways[idx] = w;
    }
    ways
}

const LETTER: [char; 4] = ['B', 'R', 'T', 'L'];

fn enumerate(pr: &Problem, ways: &[Map], idx: usize, s: &St, word: &mut Vec<char>, out: &mut dyn FnMut(&str)) {
    let n = pr.n;
    if idx == n * n {
        out(&word.iter().collect::<String>());
        return;
    }
    let mut buf = Vec::new();
    step(pr, idx, s, &mut buf);
    for (e, ns) in buf {
        if ways[idx + 1].get(&encode(&ns)).copied().unwrap_or(0) > 0 {
            word.push(LETTER[e as usize]);
            enumerate(pr, ways, idx + 1, &ns, word, out);
            word.pop();
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(|s| s.as_str()).unwrap_or("count");
    let pr = read_problem(&args[2]);
    let n = pr.n;
    match mode {
        "count" => {
            let t = std::time::Instant::now();
            let (total, _) = forward(&pr, false);
            println!("{}", total);
            eprintln!("n={} time {:.1}s", n, t.elapsed().as_secs_f64());
        }
        "stats" => {
            // per layer: reachable states, and how many of them can still complete
            let (total, layers) = forward(&pr, true);
            let ways = backward(&pr, &layers);
            let (mut sf, mut sl) = (0usize, 0usize);
            for i in 0..=n * n {
                sf += layers[i].len();
                sl += ways[i].len();
                if i % n == 0 || i == n * n / 2 {
                    println!("cell {:4}  reachable {:9}  live {:9}  ({:.1}%)", i, layers[i].len(), ways[i].len(),
                             100.0 * ways[i].len() as f64 / layers[i].len().max(1) as f64);
                }
            }
            println!("total {}  all layers: reachable {} live {} ({:.1}%)", total, sf, sl, 100.0 * sl as f64 / sf as f64);
        }
        "enum" | "sample" => {
            let (total, layers) = forward(&pr, true);
            eprintln!("total {}", total);
            let ways = backward(&pr, &layers);
            let s0 = initial(n);
            if mode == "enum" {
                let mut word = Vec::new();
                enumerate(&pr, &ways, 0, &s0, &mut word, &mut |w| println!("{}", w));
            } else {
                let k: usize = args[3].parse().unwrap();
                let mut seed: u64 = args[4].parse().unwrap();
                let mut rnd = move || {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    seed
                };
                for _ in 0..k {
                    let mut s = s0.clone();
                    let mut word = String::new();
                    for idx in 0..n * n {
                        let mut buf = Vec::new();
                        step(&pr, idx, &s, &mut buf);
                        let tot: u128 = buf.iter().map(|(_, ns)| ways[idx + 1].get(&encode(ns)).copied().unwrap_or(0)).sum();
                        let r = ((rnd() as u128) << 64 | rnd() as u128) % tot;
                        let mut acc = 0u128;
                        for (e, ns) in buf {
                            acc += ways[idx + 1].get(&encode(&ns)).copied().unwrap_or(0);
                            if r < acc {
                                word.push(LETTER[e as usize]);
                                s = ns;
                                break;
                            }
                        }
                    }
                    println!("{}", word);
                }
            }
        }
        _ => panic!("mode: count | enum | sample"),
    }
}
