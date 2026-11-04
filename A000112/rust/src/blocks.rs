//! Labeled marked-poset statistics by a three-block split, for the large Burnside quotients.
//!
//! K(N)  = #labeled posets on [N] with 1 ⋖ 2,
//! V(N)  = # with 1 ⋖ 3 and 2 ⋖ 3,            (Λ = V by duality)
//! Ch(N) = # with 1 ⋖ 2 ⋖ 3,
//! Lo(N) = # with 1 ⋖ 2 ⋖ 3 and 2 the only element strictly between 1 and 3.
//!
//! Split at the hub h (2, or 3 for V): P1 = down(h), P3 = up(h), P2 = the rest. P1 and P1 ∪ P2 are
//! ideals, so the order is (P1, P2, P3, φ: P2 → J(P1), ψ: P3 → J(P2)) with all of P1 below h and
//! P3. The marks become: 1 (and 2 for V) maximal in P1, 3 minimal in P3 for Ch and Lo. With
//! Φ_m,κ(Y) = Σ_{X labeled on [m], marks κ} Hom(X, J(Y)) (note Hom(X, J(Y)) = Hom(Y, J(X))),
//!
//!   stat(N) = Σ_{b1+b2+s=N-1} multinomial · Σ_{Y on b2} (b2!/|Aut Y|) Φ_b1,κ1(Y) Φ_s,κ3(Y).
//!
//! Lo also needs 1 ∥ b for every b ∈ P2 between 1 and 3: with U = {b : 1 ∈ φ(b)} (a filter) and
//! ψ(3) ⊆ Y∖U it becomes Σ_U Σ_{V ideal, V∩U=∅} Φ_{b1-1}(Y^{+U}) Φ_{s-1}(Y^{V+}), where Y^{+U}
//! adds a minimal point below exactly U and Y^{V+} a maximal point above exactly V.
//! (All checked by brute force for N ≤ 6 in ../README.md.)

use super::*;
use std::io::BufRead;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Base,
    Max1,
    Max2,
    Min1,
}
const KINDS: [Kind; 4] = [Kind::Base, Kind::Max1, Kind::Max2, Kind::Min1];

fn kind_idx(k: Kind) -> usize {
    match k {
        Kind::Base => 0,
        Kind::Max1 => 1,
        Kind::Max2 => 2,
        Kind::Min1 => 3,
    }
}

fn up_masks(p: &Poset) -> [u16; MAXN] {
    let mut up = [0u16; MAXN];
    for y in 0..p.n {
        for x in 0..p.n {
            if p.down[y] >> x & 1 == 1 {
                up[x] |= 1 << y;
            }
        }
    }
    up
}

fn n_max(p: &Poset) -> u32 {
    let up = up_masks(p);
    (0..p.n).filter(|&x| up[x] == 0).count() as u32
}
fn n_min(p: &Poset) -> u32 {
    (0..p.n).filter(|&x| p.down[x] == 0).count() as u32
}

/// Labeled weight of an unlabeled poset on m points with |Aut| = aut: labelings of [m] with the
/// marked labels on maximal (Max1: one, Max2: two) or minimal (Min1) elements.
fn wt(m: usize, aut: u128, mu: u32, nu: u32, k: Kind) -> u128 {
    let f = fact(m);
    let (mu, nu, m) = (mu as u128, nu as u128, m as u128);
    match k {
        Kind::Base => f / aut,
        Kind::Max1 => if m == 0 { 0 } else { f * mu / (m * aut) },
        Kind::Max2 => if m < 2 { 0 } else { f * mu * (mu.saturating_sub(1)) / (m * (m - 1) * aut) },
        Kind::Min1 => if m == 0 { 0 } else { f * nu / (m * aut) },
    }
}

// ---------------------------------------------------------------- Hom DP (from ../A006455)

struct Lattice {
    masks: Vec<u16>,
    idx: Vec<u32>,
    up: Vec<Vec<u32>>,
}

fn ideal_list(p: &Poset) -> Vec<u16> {
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

impl Lattice {
    fn new(p: &Poset) -> Self {
        let masks = ideal_list(p);
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

struct Step {
    feeds: Vec<usize>,
    new_classes: Vec<(Vec<usize>, bool)>,
}
struct Plan {
    steps: Vec<Step>,
}

impl Plan {
    fn new(s: &Poset) -> Self {
        let k = s.n;
        let cov = &s.cov[..k];
        let fut = |i: usize, j: usize| -> u16 { (j + 1..k).filter(|&t| cov[t] >> i & 1 == 1).fold(0, |m, t| m | 1 << t) };
        let mut classes: Vec<u16> = Vec::new();
        let mut steps = Vec::with_capacity(k);
        for j in 0..k {
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
            steps.push(Step { feeds, new_classes: next.into_iter().map(|e| (e.1, e.2)).collect() });
        }
        assert!(classes.is_empty());
        Plan { steps }
    }
}

/// Hom(S, L) for the plan of S (u128 counters).
fn hom_dp(plan: &Plan, lat: &Lattice, buf: &mut [Vec<u128>; 2]) -> u128 {
    let l = lat.len();
    let mut nold = 0usize;
    buf[0].clear();
    buf[0].push(1);
    let mut vals = [0usize; MAXN];
    for st in &plan.steps {
        let nnew = st.new_classes.len();
        let (a, b) = buf.split_at_mut(1);
        let (old, new) = (&a[0], &mut b[0]);
        new.clear();
        new.resize(l.pow(nnew as u32), 0);
        let j_class = st.new_classes.iter().position(|c| c.1);
        for (si, &c) in old.iter().enumerate() {
            if c == 0 {
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
            if j_class.is_none() {
                new[base] += c * ups.len() as u128;
            } else {
                for &u in ups {
                    let m = lat.masks[u as usize] | jrest;
                    new[base + lat.idx[m as usize] as usize * jmul] += c;
                }
            }
        }
        buf.swap(0, 1);
        nold = nnew;
    }
    buf[0][0]
}

/// Hom(X, J(Y)) = Hom(Y, J(X)), computed with the lattice of the smaller side.
fn hom_pair(x: &Poset, y: &Poset, buf: &mut [Vec<u128>; 2]) -> u128 {
    if x.n == 0 || y.n == 0 {
        return 1;
    }
    let (big, small) = if x.n >= y.n { (x, y) } else { (y, x) };
    hom_dp(&Plan::new(big), &Lattice::new(small), buf)
}

// ---------------------------------------------------------------- zeta functionals on J(Q)

/// Ideals of q with index map, and down/up zeta transforms over the ideal lattice.
pub struct IdealLat {
    masks: Vec<u16>,
    idx: Vec<u32>,
    up: [u16; MAXN],
    n: usize,
}

impl IdealLat {
    pub fn new(q: &Poset) -> Self {
        let masks = ideal_list(q);
        let mut idx = vec![u32::MAX; 1 << q.n];
        for (i, &m) in masks.iter().enumerate() {
            idx[m as usize] = i as u32;
        }
        IdealLat { masks, idx, up: up_masks(q), n: q.n }
    }
    pub fn len(&self) -> usize {
        self.masks.len()
    }
    /// F(I) = Σ_{J ⊆ I ideal} f(J). SOS-style: the last element processed is the first removed on
    /// the path from I down to J, so process labels upward (removal then goes top-down, through
    /// ideals only; labels are a linear extension).
    pub fn zeta_down(&self, f: &mut [u128]) {
        for x in 0..self.n {
            for (i, &m) in self.masks.iter().enumerate() {
                if m >> x & 1 == 1 && m & self.up[x] == 0 {
                    let j = self.idx[(m ^ 1 << x) as usize] as usize;
                    f[i] += f[j];
                }
            }
        }
    }
    /// F(I) = Σ_{J ⊇ I ideal} f(J). Labels processed downward, so additions go bottom-up.
    pub fn zeta_up(&self, f: &mut [u128], down: &[u16]) {
        for x in (0..self.n).rev() {
            for (i, &m) in self.masks.iter().enumerate() {
                if m >> x & 1 == 0 && m & down[x] == down[x] {
                    let j = self.idx[(m | 1 << x) as usize] as usize;
                    f[i] += f[j];
                }
            }
        }
    }
}

/// Hom(T, J(Q)) for a small forest-Hasse poset T (all posets on ≤ 3 points), by a tree DP with
/// zeta transforms: value(node, I) = product over Hasse neighbours c of Σ over I_c ⊆ I (c below)
/// or I_c ⊇ I (c above) of value(c, I_c).
pub fn hom_forest(t: &Poset, q: &Poset, lat: &IdealLat) -> u128 {
    let k = t.n;
    let l = lat.len();
    let mut adj = vec![vec![]; k]; // (neighbour, neighbour below?)
    for y in 0..k {
        for x in 0..k {
            if t.cov[y] >> x & 1 == 1 {
                adj[y].push((x, true));
                adj[x].push((y, false));
            }
        }
    }
    let mut seen = vec![false; k];
    let mut total: u128 = 1;
    fn val(v: usize, parent: usize, adj: &[Vec<(usize, bool)>], q: &Poset, lat: &IdealLat, seen: &mut [bool]) -> Vec<u128> {
        seen[v] = true;
        let mut f = vec![1u128; lat.len()];
        for &(c, below) in &adj[v] {
            if c == parent {
                continue;
            }
            let mut g = val(c, v, adj, q, lat, seen);
            if below {
                lat.zeta_down(&mut g);
            } else {
                lat.zeta_up(&mut g, &q.down);
            }
            for i in 0..f.len() {
                f[i] *= g[i];
            }
        }
        f
    }
    for r in 0..k {
        if !seen[r] {
            let f = val(r, usize::MAX, &adj, q, lat, &mut seen);
            total *= f.iter().sum::<u128>();
        }
    }
    let _ = l;
    total
}

// ---------------------------------------------------------------- families and aggregates

pub struct Fam {
    pub ps: Vec<Poset>,
    pub w: Vec<[u128; 4]>, // labeled weights per kind
}

fn fam(m: usize) -> Fam {
    let ps = load(m);
    let w = ps
        .iter()
        .map(|p| {
            let a = automorphisms(p);
            let (mu, nu) = (n_max(p), n_min(p));
            [wt(m, a, mu, nu, Kind::Base), wt(m, a, mu, nu, Kind::Max1), wt(m, a, mu, nu, Kind::Max2), wt(m, a, mu, nu, Kind::Min1)]
        })
        .collect();
    Fam { ps, w }
}

/// Are p and q isomorphic (tiny posets only)?
fn iso(p: &Poset, q: &Poset) -> bool {
    if p.n != q.n {
        return false;
    }
    let n = p.n;
    let mut perm: Vec<usize> = (0..n).collect();
    fn rec(i: usize, n: usize, p: &Poset, q: &Poset, perm: &mut Vec<usize>, used: u16) -> bool {
        if i == n {
            return true;
        }
        for v in 0..n {
            if used >> v & 1 == 1 {
                continue;
            }
            perm[i] = v;
            if (0..i).all(|j| (p.down[i] >> j & 1) == (q.down[v] >> perm[j] & 1) && (p.down[j] >> i & 1) == (q.down[perm[j]] >> v & 1))
                && rec(i + 1, n, p, q, perm, used | 1 << v)
            {
                return true;
            }
        }
        false
    }
    rec(0, n, p, q, &mut perm, 0)
}

/// The tiny posets (≤ 3 points), indexed as tiny[k][i].
pub struct Tiny {
    pub fams: Vec<Fam>,
}
impl Tiny {
    fn new() -> Self {
        Tiny { fams: (0..=3).map(fam).collect() }
    }
    fn index(&self, z: &Poset) -> usize {
        self.fams[z.n].ps.iter().position(|t| iso(t, z)).expect("tiny iso")
    }
}

/// Aggregates of a streamed family (m ≥ 11): ph[kind][k][i] = Σ_Q wt_kind(Q) Hom(tiny[k][i], J(Q))
/// and prod[a][ka][b][kb] = Σ_Q w0(Q) Φ_a,ka(Q) Φ_b,kb(Q) for tiny a, b, plus G[k] = Σ_Q w0 d^k.
#[derive(Default, Clone)]
pub struct Agg {
    pub m: usize,
    pub count: u64,
    pub ph: Vec<Vec<Vec<u128>>>,
    pub prod: Vec<u128>,
    pub g: Vec<u128>,
}

fn prod_index(a: usize, ka: usize, b: usize, kb: usize) -> usize {
    ((a * 4 + ka) * 4 + b) * 4 + kb
}

impl Agg {
    fn empty(m: usize, tiny: &Tiny) -> Self {
        Agg {
            m,
            count: 0,
            ph: (0..4).map(|_| (0..=3).map(|k| vec![0; tiny.fams[k].ps.len()]).collect()).collect(),
            prod: vec![0; 4 * 4 * 4 * 4],
            g: vec![0; 6],
        }
    }
    fn to_text(&self) -> String {
        let mut s = format!("m {}\ncount {}\n", self.m, self.count);
        for kd in 0..4 {
            for k in 0..=3 {
                let v: Vec<String> = self.ph[kd][k].iter().map(|x| x.to_string()).collect();
                s += &format!("ph {kd} {k} {}\n", v.join(" "));
            }
        }
        let v: Vec<String> = self.prod.iter().map(|x| x.to_string()).collect();
        s += &format!("prod {}\n", v.join(" "));
        let v: Vec<String> = self.g.iter().map(|x| x.to_string()).collect();
        s += &format!("g {}\n", v.join(" "));
        s
    }
    fn add_text(&mut self, text: &str) {
        for line in text.lines() {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.is_empty() {
                continue;
            }
            match f[0] {
                "m" => assert_eq!(f[1].parse::<usize>().unwrap(), self.m),
                "count" => self.count += f[1].parse::<u64>().unwrap(),
                "ph" => {
                    let (kd, k): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
                    for (i, x) in f[3..].iter().enumerate() {
                        self.ph[kd][k][i] += x.parse::<u128>().unwrap();
                    }
                }
                "prod" => {
                    for (i, x) in f[1..].iter().enumerate() {
                        self.prod[i] += x.parse::<u128>().unwrap();
                    }
                }
                "g" => {
                    for (i, x) in f[1..].iter().enumerate() {
                        self.g[i] += x.parse::<u128>().unwrap();
                    }
                }
                _ => panic!("bad agg line {line}"),
            }
        }
    }
}

/// Per-poset tiny data: h[k][i] = Hom(tiny[k][i], J(q)).
fn tiny_homs(q: &Poset, tiny: &Tiny) -> Vec<Vec<u128>> {
    let lat = IdealLat::new(q);
    (0..=3).map(|k| tiny.fams[k].ps.iter().map(|t| hom_forest(t, q, &lat)).collect()).collect()
}

/// Φ_a,kind(Q) for tiny a from the tiny homs of Q.
fn phi_tiny(h: &[Vec<u128>], tiny: &Tiny, a: usize, kd: usize) -> u128 {
    tiny.fams[a].w.iter().zip(&h[a]).map(|(w, x)| w[kd] * x).sum()
}

fn agg_add(agg: &mut Agg, q: &Poset, tiny: &Tiny) {
    let m = q.n;
    let aut = automorphisms(q);
    let (mu, nu) = (n_max(q), n_min(q));
    let w: Vec<u128> = KINDS.iter().map(|&k| wt(m, aut, mu, nu, k)).collect();
    let h = tiny_homs(q, tiny);
    for kd in 0..4 {
        for k in 0..=3 {
            for i in 0..h[k].len() {
                agg.ph[kd][k][i] += w[kd] * h[k][i];
            }
        }
    }
    let mut phi = [[0u128; 4]; 4];
    for a in 0..=3 {
        for kd in 0..4 {
            phi[a][kd] = phi_tiny(&h, tiny, a, kd);
        }
    }
    for a in 0..=3 {
        for ka in 0..4 {
            for b in 0..=3 {
                for kb in 0..4 {
                    if a + b <= 3 {
                        agg.prod[prod_index(a, ka, b, kb)] += w[0] * phi[a][ka] * phi[b][kb];
                    }
                }
            }
        }
    }
    let d = h[1][0];
    let mut p = 1u128;
    for k in 0..6 {
        agg.g[k] += w[0] * p;
        p = p.saturating_mul(d);
    }
    agg.count += 1;
}

/// `blocks agg M`: aggregate the M-point posets on stdin (genposetg M t); prints the aggregate.
pub fn cmd_agg(m: usize) {
    disable_power_throttling();
    let t0 = Instant::now();
    let tiny = Tiny::new();
    let mut agg = Agg::empty(m, &tiny);
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        if line.starts_with('&') {
            let q = parse_d6(&line);
            assert_eq!(q.n, m);
            agg_add(&mut agg, &q, &tiny);
        }
    }
    eprintln!("agg {m}: {} posets, {:.1}s", agg.count, t0.elapsed().as_secs_f64());
    print!("{}", agg.to_text());
}

// ---------------------------------------------------------------- the engine

pub struct Engine {
    pub fams: Vec<Fam>, // m ≤ 10
    pub tiny: Tiny,
    pub aggs: BTreeMap<usize, Agg>, // m = 11, 12
    pub g: BTreeMap<(usize, usize), u128>,
    buf: [Vec<u128>; 2],
}

impl Engine {
    pub fn new(maxmem: usize) -> Self {
        let t0 = Instant::now();
        let fams: Vec<Fam> = (0..=maxmem).map(fam).collect();
        let tiny = Tiny::new();
        let mut aggs = BTreeMap::new();
        for m in [11usize, 12] {
            let path = format!("{}/../aggs/agg{m}.txt", env!("CARGO_MANIFEST_DIR"));
            if let Ok(text) = std::fs::read_to_string(&path) {
                let mut a = Agg::empty(m, &tiny);
                a.add_text(&text);
                let lab: u128 = a.ph[0][0][0];
                assert_eq!(lab, A001035[m], "agg{m}: labeled count mismatch (incomplete?)");
                aggs.insert(m, a);
            }
        }
        eprintln!("engine: families ≤ {maxmem}, aggregates {:?} ({:.1}s)", aggs.keys().collect::<Vec<_>>(), t0.elapsed().as_secs_f64());
        let mut e = Engine { fams, tiny, aggs, g: BTreeMap::new(), buf: [Vec::new(), Vec::new()] };
        e.fill_g();
        e
    }

    /// G(m, k) = Σ_{labeled Q on m} d(Q)^k for the m available, and G(13,1) from P(14).
    fn fill_g(&mut self) {
        for m in 0..self.fams.len() {
            for k in 0..=14 - m {
                let s: u128 = self.fams[m]
                    .ps
                    .iter()
                    .zip(&self.fams[m].w)
                    .map(|(p, w)| w[0] * (ideal_list(p).len() as u128).pow(k as u32))
                    .sum();
                self.g.insert((m, k), s);
            }
        }
        for (&m, a) in &self.aggs {
            for k in 0..6 {
                self.g.insert((m, k), a.g[k]);
            }
        }
        // P(N) = Σ_{k≥1} (-1)^(k+1) C(N,k) G(N-k,k)  (peeling all maximal elements)
        for nn in [12usize, 13, 14] {
            let target = (nn - 1, 1);
            if self.g.contains_key(&target) {
                continue;
            }
            let mut acc: i128 = A001035[nn] as i128;
            let mut ok = true;
            for k in 2..=nn {
                match self.g.get(&(nn - k, k)) {
                    Some(&v) => {
                        let c = binom(nn, k) as i128 * v as i128;
                        if k % 2 == 0 { acc += c } else { acc -= c }
                    }
                    None => ok = false,
                }
            }
            if ok {
                assert!(acc % nn as i128 == 0);
                self.g.insert(target, (acc / nn as i128) as u128);
            }
        }
    }

    pub fn check_g(&self) {
        // the identity also holds where both sides are known
        for nn in 1..self.fams.len() {
            let mut acc: i128 = 0;
            for k in 1..=nn {
                let c = binom(nn, k) as i128 * self.g[&(nn - k, k)] as i128;
                if k % 2 == 1 { acc += c } else { acc -= c }
            }
            assert_eq!(acc as u128, A001035[nn], "peeling identity at {nn}");
        }
    }

    fn gv(&self, m: usize, k: usize) -> u128 {
        *self.g.get(&(m, k)).unwrap_or_else(|| panic!("G({m},{k}) unavailable"))
    }

    /// Φ_m,kind(Z) = Σ_{X labeled on m, marks} Hom(X, J(Z)) for one poset Z.
    pub fn phi(&mut self, m: usize, kd: Kind, z: &Poset) -> u128 {
        let ki = kind_idx(kd);
        if z.n == 0 {
            // labeled counts with marks
            return match kd {
                Kind::Base => A001035[m],
                Kind::Max1 | Kind::Min1 => if m == 0 { 0 } else { self.gv(m - 1, 1) },
                Kind::Max2 => if m < 2 { 0 } else { self.gv(m - 2, 2) },
            };
        }
        if m < self.fams.len() {
            let mut buf = std::mem::take(&mut self.buf);
            let mut s = 0u128;
            for i in 0..self.fams[m].ps.len() {
                let w = self.fams[m].w[i][ki];
                if w != 0 {
                    s += w * hom_pair(&self.fams[m].ps[i], z, &mut buf);
                }
            }
            self.buf = buf;
            return s;
        }
        if let Some(a) = self.aggs.get(&m) {
            assert!(z.n <= 3, "streamed family {m} needs a tiny partner, got {}", z.n);
            let i = self.tiny.index(z);
            return a.ph[ki][z.n][i];
        }
        // m = 13 with one point: Σ_{X on 13, 1 max} d(X) = G(12,2) + H3(12)
        if m == 13 && kd == Kind::Max1 && z.n == 1 {
            let a = &self.aggs[&12];
            let chain2 = self.tiny.index(&chain(2));
            return self.gv(12, 2) + a.ph[0][2][chain2];
        }
        panic!("Φ_{m},{kd:?} on {} points unavailable", z.n)
    }

    /// Vector Φ_m,kind(Y) over the family of size b2 (in memory).
    fn phi_vec(&mut self, m: usize, kd: Kind, b2: usize) -> Vec<u128> {
        let ki = kind_idx(kd);
        let ny = self.fams[b2].ps.len();
        if m == 0 {
            return vec![if kd == Kind::Base { 1 } else { 0 }; ny];
        }
        if m < self.fams.len() {
            let mut buf = std::mem::take(&mut self.buf);
            let mut out = vec![0u128; ny];
            let (xs, ys) = (&self.fams[m], &self.fams[b2]);
            if b2 == 0 {
                out[0] = xs.w.iter().map(|w| w[ki]).sum();
            } else {
                // Hom(X, J(Y)) = Hom(Y, J(X)): the dense DP state of one orientation is
                // |J(other)|^width(plan); flip a pair when its state would be large and the
                // other way is cheaper (wide 9- and 10-point X against 5- to 7-point Y).
                let side = |f: &Fam| -> (Vec<Plan>, Vec<u32>, Vec<f64>) {
                    let plans: Vec<Plan> = f.ps.iter().map(Plan::new).collect();
                    let w = plans.iter().map(|p| p.steps.iter().map(|s| s.new_classes.len()).max().unwrap_or(0) as u32).collect();
                    let l = f.ps.iter().map(|p| ideal_list(p).len() as f64).collect();
                    (plans, w, l)
                };
                let (xp, xw, xl) = side(xs);
                let (yp, yw, yl) = side(ys);
                let ylat: Vec<Lattice> = if m >= b2 { ys.ps.iter().map(Lattice::new).collect() } else { vec![] };
                let xlat: Vec<Lattice> = if m < b2 { xs.ps.iter().map(Lattice::new).collect() } else { vec![] };
                const BIG: f64 = (1u64 << 20) as f64;
                for (i, w) in xs.w.iter().enumerate() {
                    if w[ki] == 0 {
                        continue;
                    }
                    if xs.ps[i].n == 0 {
                        for o in out.iter_mut() {
                            *o += w[ki];
                        }
                        continue;
                    }
                    let mut own: Option<Lattice> = None;
                    for j in 0..ny {
                        // default: plan of the larger side, lattice of the smaller
                        let fwd = m >= b2;
                        let cost_xy = yl[j].powi(xw[i] as i32); // plan X, lattice J(Y)
                        let cost_yx = xl[i].powi(yw[j] as i32); // plan Y, lattice J(X)
                        let use_xy = if fwd { !(cost_xy > BIG && cost_yx < cost_xy) } else { cost_yx > BIG && cost_xy < cost_yx };
                        let h = if use_xy {
                            let lat = if fwd { &ylat[j] } else { &*own.get_or_insert_with(|| Lattice::new(&ys.ps[j])) };
                            hom_dp(&xp[i], lat, &mut buf)
                        } else {
                            let lat = if fwd { &*own.get_or_insert_with(|| Lattice::new(&xs.ps[i])) } else { &xlat[i] };
                            hom_dp(&yp[j], lat, &mut buf)
                        };
                        if !fwd {
                            own = None;
                        }
                        out[j] += w[ki] * h;
                    }
                }
            }
            self.buf = buf;
            return out;
        }
        let ys: Vec<Poset> = self.fams[b2].ps.clone();
        ys.iter().map(|y| self.phi(m, kd, y)).collect()
    }

    /// K, Ch, V: Σ over block sizes of multinomial · Σ_Y w0(Y) Φ_b1,κ1(Y) Φ_s,κ3(Y).
    pub fn stat_direct(&mut self, n: usize, k1: Kind, k3: Kind) -> u128 {
        let (m1, m3) = (if k1 == Kind::Max2 { 2 } else { 1 }, if k3 == Kind::Min1 { 1 } else { 0 });
        let mut total = 0u128;
        for b1 in m1..n {
            for s in m3..n - b1 {
                let b2 = n - 1 - b1 - s;
                let mult = fact(n - 1 - m1 - m3) / (fact(b1 - m1) * fact(b2) * fact(s - m3));
                let inner: u128 = if b2 < self.fams.len() {
                    let a = self.phi_vec(b1, k1, b2);
                    let c = self.phi_vec(s, k3, b2);
                    (0..a.len()).map(|j| self.fams[b2].w[j][0] * a[j] * c[j]).sum()
                } else if let Some(ag) = self.aggs.get(&b2) {
                    assert!(b1 + s <= 3);
                    ag.prod[prod_index(b1, kind_idx(k1), s, kind_idx(k3))]
                } else if b2 == 13 && b1 == 1 && s == 0 && k1 == Kind::Max1 && k3 == Kind::Base {
                    self.gv(13, 1)
                } else {
                    panic!("no data for Y on {b2} points (b1={b1}, s={s})")
                };
                total += mult * inner;
            }
        }
        total
    }

    /// Lo(n) by the U/V decomposition.
    pub fn stat_lo(&mut self, n: usize) -> u128 {
        let mut total = 0u128;
        for b1 in 1..n {
            for s in 1..n - b1 {
                let b2 = n - 1 - b1 - s;
                let mult = fact(n - 3) / (fact(b1 - 1) * fact(b2) * fact(s - 1));
                let inner: u128 = if b2 <= 8 && b2 < self.fams.len() {
                    // generic: explicit filters U and ideals V
                    let ys = self.fams[b2].ps.clone();
                    let ws: Vec<u128> = self.fams[b2].w.iter().map(|w| w[0]).collect();
                    let mut acc = 0u128;
                    for (y, w0) in ys.iter().zip(ws) {
                        let lat = IdealLat::new(y);
                        let full: u16 = ((1u32 << y.n) - 1) as u16;
                        // a(I) = Φ_{b1-1}(Y^{+U}) with U = complement of I; b(V) = Φ_{s-1}(Y^{V+})
                        let mut a: Vec<u128> = lat.masks.iter().map(|&i| self.phi(b1 - 1, Kind::Base, &add_below(y, full & !i))).collect();
                        let bv: Vec<u128> = lat.masks.iter().map(|&v| self.phi(s - 1, Kind::Base, &add_above(y, v))).collect();
                        lat.zeta_up(&mut a, &y.down); // Σ_{I ⊇ V} a(I) = Σ over filters U ⊆ complement(V)
                        acc += w0 * (0..a.len()).map(|i| a[i] * bv[i]).sum::<u128>();
                    }
                    acc
                } else if b2 == 9 || b2 == 10 {
                    assert!(b1 - 1 <= 2 && s - 1 <= 2);
                    let (ys, ws): (Vec<Poset>, Vec<u128>) = if b2 < self.fams.len() {
                        (self.fams[b2].ps.clone(), self.fams[b2].w.iter().map(|w| w[0]).collect())
                    } else {
                        panic!("Lo needs posets on {b2} points in memory")
                    };
                    let mut acc = 0u128;
                    for (y, w0) in ys.iter().zip(ws) {
                        acc += w0 * lo_closed(y, b1 - 1, s - 1);
                    }
                    acc
                } else if b2 == 11 && b1 == 1 && s == 1 {
                    let chain2 = self.tiny.index(&chain(2));
                    self.aggs[&11].ph[0][2][chain2]
                } else {
                    panic!("Lo: no data for Y on {b2} points (b1={b1}, s={s})")
                };
                total += mult * inner;
            }
        }
        total
    }
}

/// Closed form of Σ_U Σ_{V∩U=∅} Φ_a(Y^{+U}) Φ_b(Y^{V+}) for a, b ≤ 2 (labeled X on ≤ 2 points:
/// Φ_0 = 1, Φ_1(Z) = d(Z), Φ_2(Z) = d(Z)^2 + 2 h3(Z), h3 = nested ideal pairs).
fn lo_closed(y: &Poset, a: usize, b: usize) -> u128 {
    let lat = IdealLat::new(y);
    let l = lat.len();
    let full: u16 = ((1u32 << y.n) - 1) as u16;
    let d = l as u128;
    let mut down = vec![1u128; l];
    lat.zeta_down(&mut down); // #ideals ⊆ I
    let mut up = vec![1u128; l];
    lat.zeta_up(&mut up, &y.down); // #ideals ⊇ I
    let h3: u128 = down.iter().sum();
    let mut zd_up = up.clone();
    lat.zeta_down(&mut zd_up); // Σ_{J ⊆ I} up(J)
    let mut zd_down = down.clone();
    lat.zeta_down(&mut zd_down); // Σ_{J ⊆ I} down(J)
    let mut zu_down = down.clone();
    lat.zeta_up(&mut zu_down, &y.down); // Σ_{J ⊇ I} down(J)
    let mut zu_up = up.clone();
    lat.zeta_up(&mut zu_up, &y.down); // Σ_{J ⊇ I} up(J)
    // For U = complement of ideal I (Uc = I): Z = Y^{+U}: d(Z) = d + down(I);
    // h3(Z) = h3 + Σ_{J ⊆ I} up(J) + Σ_{J ⊆ I} down(J).
    let phi_u = |i: usize| -> u128 {
        let dz = d + down[i];
        match a {
            0 => 1,
            1 => dz,
            2 => dz * dz + 2 * (h3 + zd_up[i] + zd_down[i]),
            _ => unreachable!(),
        }
    };
    // For V = ideal i: Z = Y^{V+}: d(Z) = d + up(V); h3(Z) = h3 + Σ_{J⊇V} down(J) + Σ_{J⊇V} up(J).
    let phi_v = |i: usize| -> u128 {
        let dz = d + up[i];
        match b {
            0 => 1,
            1 => dz,
            2 => dz * dz + 2 * (h3 + zu_down[i] + zu_up[i]),
            _ => unreachable!(),
        }
    };
    let mut av: Vec<u128> = (0..l).map(phi_u).collect();
    lat.zeta_up(&mut av, &y.down);
    let _ = full;
    (0..l).map(|i| av[i] * phi_v(i)).sum()
}

fn chain(k: usize) -> Poset {
    let mut p = Poset { n: k, down: [0; MAXN], cov: [0; MAXN] };
    for j in 1..k {
        p.cov[j] = 1 << (j - 1);
        p.down[j] = (1 << j) - 1;
    }
    p
}

fn closure(p: &mut Poset) {
    // recompute cov from down (down must be transitively closed)
    for y in 0..p.n {
        let mut c = p.down[y];
        for x in 0..p.n {
            if p.down[y] >> x & 1 == 1 {
                c &= !p.down[x];
            }
        }
        p.cov[y] = c;
    }
}

/// Y plus a new minimal point below exactly the filter U, relabeled so that labels stay a linear
/// extension (the new point gets label 0).
fn add_below(y: &Poset, u: u16) -> Poset {
    let n = y.n;
    let mut p = Poset { n: n + 1, down: [0; MAXN], cov: [0; MAXN] };
    for j in 0..n {
        p.down[j + 1] = (y.down[j] << 1) | if u >> j & 1 == 1 { 1 } else { 0 };
    }
    closure(&mut p);
    p
}

/// Y plus a new maximal point above exactly the ideal V (label n).
fn add_above(y: &Poset, v: u16) -> Poset {
    let n = y.n;
    let mut p = Poset { n: n + 1, down: [0; MAXN], cov: [0; MAXN] };
    p.down[..n].copy_from_slice(&y.down[..n]);
    p.down[n] = v;
    closure(&mut p);
    p
}

pub fn binom(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    let mut r = 1u128;
    for i in 0..k {
        r = r * (n - i) as u128 / (i + 1) as u128;
    }
    r
}

/// `blocks N`: K, V, Ch, Lo for 2..=N (as far as the data reaches).
pub fn cmd_blocks(nmax: usize, maxmem: usize) {
    disable_power_throttling();
    let t0 = Instant::now();
    let mut e = Engine::new(maxmem);
    e.check_g();
    // zeta-based tiny homs agree with the DP on every poset up to 8 points
    let mut buf = [Vec::new(), Vec::new()];
    for m in 0..=8.min(maxmem) {
        for q in &e.fams[m].ps {
            let h = tiny_homs(q, &e.tiny);
            for k in 0..=3 {
                for (i, t) in e.tiny.fams[k].ps.iter().enumerate() {
                    assert_eq!(h[k][i], hom_pair(t, q, &mut buf), "tiny hom m={m}");
                }
            }
        }
    }
    eprintln!("tiny homs checked (m <= 8)");
    for n in 2..=nmax {
        let k = e.stat_direct(n, Kind::Max1, Kind::Base);
        let (v, ch, lo) = if n >= 3 {
            (e.stat_direct(n, Kind::Max2, Kind::Base), e.stat_direct(n, Kind::Max1, Kind::Min1), e.stat_lo(n))
        } else {
            (0, 0, 0)
        };
        println!("N={n} K={k} V={v} Ch={ch} Lo={lo}  ({:.1}s)", t0.elapsed().as_secs_f64());
    }
}

/// `blocks1 N WHICH`: one statistic (K, V, Ch or Lo) at one N.
pub fn cmd_blocks1(n: usize, which: &str, maxmem: usize) {
    disable_power_throttling();
    let t0 = Instant::now();
    let mut e = Engine::new(maxmem);
    let v = match which {
        "K" => e.stat_direct(n, Kind::Max1, Kind::Base),
        "V" => e.stat_direct(n, Kind::Max2, Kind::Base),
        "Ch" => e.stat_direct(n, Kind::Max1, Kind::Min1),
        "Lo" => e.stat_lo(n),
        _ => panic!("which = K, V, Ch or Lo"),
    };
    println!("{which}({n}) = {v}  ({:.1}s)", t0.elapsed().as_secs_f64());
}
