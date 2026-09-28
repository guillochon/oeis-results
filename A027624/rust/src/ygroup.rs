//! "Union grouping": a(d+2) = sum over vertex sets Y of Q_d of c(Y) * f(~Y)^2.
//!
//! In the pair sum a(d+2) = sum_{S,U} f(~(S | U))^2, group the pairs by Y = S | U. The number of
//! pairs of independent sets with union Y is
//!
//!     c(Y) = 3^(isolated vertices of Q_d[Y]) * 2^(other components of Q_d[Y]),
//!
//! since an isolated vertex can be in S, U or both, and every larger component (bipartite,
//! connected) has exactly two proper 2-colourings into S-only / U-only.
//!
//! Layout: Q_d = Q_m x C_4 with m = d - 2, i.e. four fibres of k = 2^m bits, Y = Y0|Y1|Y2|Y3.
//! Fibres 0,1 form the "lo" copy of Q_{d-1} and fibres 3,2 the "hi" copy (Q_d = Q_{d-1} x K_2).
//! Y_lo = (Y0, Y1) runs over orbit representatives of Aut(Q_{d-1}) (acting on both copies at
//! once), Y3 and Y2 over everything. For fixed (Y0, Y1, Y3) the whole slice W2 -> f(W) is one
//! subset-sum transform over k bits:
//!
//!     f(W) = sum_{t2 in I_m, t2 inside W2} G(t2),
//!     G(t2) = sum_{t0 in I_m, t0 inside W0} f_m(W1 \ (t0|t2)) * f_m(W3 \ (t0|t2)).
//!
//! Sums are exact in 256 bits (a(8) is about 2^130).

use crate::{generators, indep_sets, orbits, Group};
#[allow(unused_imports)]
use crate::FEval;
use std::time::Instant;

/// Unsigned 256-bit accumulator.
#[derive(Clone, Copy, Default, PartialEq, Debug)]
pub struct U256 {
    pub hi: u128,
    pub lo: u128,
}

impl U256 {
    pub fn add(&mut self, o: U256) {
        let (lo, carry) = self.lo.overflowing_add(o.lo);
        self.lo = lo;
        self.hi = self.hi.wrapping_add(o.hi).wrapping_add(carry as u128);
    }
    /// Full 128 x 128 -> 256-bit product.
    pub fn mul(a: u128, b: u128) -> U256 {
        let (a1, a0) = (a >> 64, a & (u64::MAX as u128));
        let (b1, b0) = (b >> 64, b & (u64::MAX as u128));
        let (p00, p01, p10, p11) = (a0 * b0, a0 * b1, a1 * b0, a1 * b1);
        let mid = (p00 >> 64) + (p01 & (u64::MAX as u128)) + (p10 & (u64::MAX as u128));
        let lo = (p00 & (u64::MAX as u128)) | (mid << 64);
        let hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
        U256 { hi, lo }
    }
    pub fn scale(self, w: u64) -> U256 {
        let mut r = U256::mul(self.lo, w as u128);
        r.hi = r.hi.wrapping_add(self.hi.wrapping_mul(w as u128));
        r
    }
    pub fn to_decimal(self) -> String {
        if self.hi == 0 {
            return self.lo.to_string();
        }
        // repeated division by 10^19 on four 64-bit limbs
        let mut limbs = [(self.lo & u64::MAX as u128) as u64, (self.lo >> 64) as u64, (self.hi & u64::MAX as u128) as u64, (self.hi >> 64) as u64];
        let mut parts = vec![];
        while limbs.iter().any(|&x| x != 0) {
            let mut rem: u128 = 0;
            for l in limbs.iter_mut().rev() {
                let cur = (rem << 64) | *l as u128;
                *l = (cur / 10_000_000_000_000_000_000) as u64;
                rem = cur % 10_000_000_000_000_000_000;
            }
            parts.push(rem as u64);
        }
        let mut s = parts.last().unwrap().to_string();
        for p in parts.iter().rev().skip(1) {
            s += &format!("{:019}", p);
        }
        s
    }
}

/// Everything that depends only on d.
pub struct YGroup {
    #[allow(dead_code)]
    pub d: usize,
    k: usize,
    km: u64,     // mask of one fibre
    full: u64,   // mask of all 4k bits
    bitmask: Vec<u64>, // bitmask[j]: vertices whose fibre position has bit j clear
    fm: Vec<u32>,      // fm[X] = # independent sets of Q_m inside X (k-bit X)
    im: Vec<u32>,      // independent sets of Q_m
    off: Vec<u32>,     // CSR over X: independent sets of Q_m inside X
    lst: Vec<u32>,
    pow3: Vec<u128>,
    pow2: Vec<u128>,
    qnb: Vec<u32>,     // qnb[p] = Q_m neighbours of fibre position p
    isoq: Vec<u32>,    // isoq[X] = positions of X with no Q_m neighbour in X
}

/// Per-slice data for the fast c(Y): everything that depends on the fixed fibres 0, 1, 3.
pub struct SliceCtx {
    nb: [[u32; 256]; 2], // extended neighbours of the fibre-2 positions, by byte
    y13: u32,            // positions occupied in fibre 1 or 3
    i1: u32,             // positions of fixed vertices isolated in the fixed part, fibre 1
    i3: u32,             // ... fibre 3
    base_iso: u32,       // isolated fixed vertices of fibre 0 (they stay isolated)
    base_c2: u32,        // fixed components of size >= 2 not touching fibre 2 positions
    touch2: Vec<u32>,    // contact masks of the other fixed components of size >= 2
}

/// Bins of the (isolated vertices, components) histogram.
pub const MAXI: usize = 65;
pub const MAXC: usize = 34; // components <= 32, plus one for the swap-symmetry weight 2

impl YGroup {
    pub fn new(d: usize) -> YGroup {
        assert!((3..=6).contains(&d));
        let m = d - 2;
        let k = 1usize << m;
        let km = if k == 64 { u64::MAX } else { (1u64 << k) - 1 };
        let full = if 4 * k == 64 { u64::MAX } else { (1u64 << (4 * k)) - 1 };
        let bitmask = (0..m)
            .map(|j| (0..4 * k).filter(|p| p >> j & 1 == 0).fold(0u64, |a, p| a | 1 << p))
            .collect();
        let im = indep_sets(m);
        let size = 1usize << k;
        let mut fm = vec![0u32; size];
        for &a in &im {
            fm[a as usize] = 1;
        }
        for b in 0..k {
            for x in 0..size {
                if x >> b & 1 == 1 {
                    fm[x] += fm[x ^ 1 << b];
                }
            }
        }
        let mut off = vec![0u32; size + 1];
        for x in 0..size {
            off[x + 1] = off[x] + im.iter().filter(|&&a| a as usize & !x == 0).count() as u32;
        }
        let mut lst = Vec::with_capacity(off[size] as usize);
        for x in 0..size {
            lst.extend(im.iter().filter(|&&a| a as usize & !x == 0));
        }
        let pow3 = (0..=64).map(|i| 3u128.pow(i)).collect();
        let pow2 = (0..=64).map(|i| 1u128 << i).collect();
        let qnb: Vec<u32> = (0..k).map(|p| (0..m).fold(0u32, |a, j| a | 1 << (p ^ 1 << j))).collect();
        let isoq = (0..size as u32)
            .map(|x| {
                let n = (0..k).filter(|p| x >> p & 1 == 1).fold(0u32, |a, p| a | qnb[p]);
                x & !n
            })
            .collect();
        YGroup { d, k, km, full, bitmask, fm, im, off, lst, pow3, pow2, qnb, isoq }
    }

    /// Precomputes the fast c(Y) data for fixed fibres 0, 1, 3.
    pub fn ctx(&self, y0: u64, y1: u64, y3: u64) -> SliceCtx {
        let k = self.k;
        let fixed = y0 | y1 << k | y3 << (3 * k);
        let iso_f = fixed & !self.nbrs(fixed);
        let km = self.km;
        let mut c = SliceCtx {
            nb: [[0; 256]; 2],
            y13: ((y1 | y3) & km) as u32,
            i1: (iso_f >> k & km) as u32,
            i3: (iso_f >> (3 * k) & km) as u32,
            base_iso: (iso_f & km).count_ones(),
            base_c2: 0,
            touch2: vec![],
        };
        let mut a: Vec<u32> = self.qnb.clone();
        let mut rem = fixed & !iso_f;
        while rem != 0 {
            let mut comp = rem & rem.wrapping_neg();
            loop {
                let nxt = (comp | self.nbrs(comp)) & fixed;
                if nxt == comp {
                    break;
                }
                comp = nxt;
            }
            rem &= !comp;
            let tau = ((comp >> k | comp >> (3 * k)) & km) as u32;
            if tau == 0 {
                c.base_c2 += 1;
            } else {
                c.touch2.push(tau);
                for p in 0..k {
                    if tau >> p & 1 == 1 {
                        a[p] |= tau;
                    }
                }
            }
        }
        for (ch, tab) in c.nb.iter_mut().enumerate() {
            for (b, e) in tab.iter_mut().enumerate() {
                for j in 0..8 {
                    let p = 8 * ch + j;
                    if b >> j & 1 == 1 && p < k {
                        *e |= a[p];
                    }
                }
            }
        }
        c
    }

    /// (isolated vertices, components of size >= 2) of Q_d[Y] for Y = fixed fibres + y2.
    #[inline]
    pub fn c_fast(&self, cx: &SliceCtx, y2: u32) -> (usize, usize) {
        let isoq = self.isoq[y2 as usize] & !cx.y13;
        let iso = cx.base_iso + (cx.i1 & !y2).count_ones() + (cx.i3 & !y2).count_ones() + isoq.count_ones();
        // components of y2 under the extended adjacency (each absorbs the fixed components it touches)
        let mut rem = y2 & !isoq;
        let mut n = cx.base_c2;
        while rem != 0 {
            let mut comp = rem & rem.wrapping_neg();
            loop {
                let ext = cx.nb[0][(comp & 255) as usize] | cx.nb[1][(comp >> 8 & 255) as usize];
                let nxt = (comp | ext) & y2;
                if nxt == comp {
                    break;
                }
                comp = nxt;
            }
            rem &= !comp;
            n += 1;
        }
        for &t in &cx.touch2 {
            n += (t & y2 == 0) as u32;
        }
        (iso as usize, n as usize)
    }

    /// Combines an (iso, comps) histogram of f^2 sums into sum c(Y) f^2.
    pub fn combine(&self, hist: &[u128]) -> U256 {
        let mut acc = U256::default();
        for i in 0..MAXI {
            for c in 0..MAXC {
                let h = hist[i * MAXC + c];
                if h != 0 {
                    acc.add(U256::mul(self.pow3[i] * self.pow2[c], h));
                }
            }
        }
        acc
    }

    /// Fast version of rep_sum: c(Y) through SliceCtx, f^2 accumulated per (iso, comps) bin.
    ///
    /// With `sym`, the lo <-> hi swap symmetry is used: the summand is invariant under swapping
    /// Y_lo and Y_hi, and popcount is invariant under Aut(Q_{d-1}), so terms with
    /// |Y_hi| < |Y_lo| are skipped, |Y_hi| = |Y_lo| count once and |Y_hi| > |Y_lo| count twice.
    /// The weighted rep sums differ from the plain ones, but their weighted total is the same.
    pub fn rep_sum_fast(&self, y_lo: u64, y3_range: std::ops::Range<u64>, sym: bool) -> U256 {
        let k = self.k;
        let (y0, y1) = (y_lo & self.km, y_lo >> k & self.km);
        let pl = y_lo.count_ones();
        let mut phi = vec![0u64; 1usize << k];
        let mut hist = vec![0u128; MAXI * MAXC]; // per bin <= 2^33 * f^2 < 2^104
        for y3 in y3_range {
            if sym && y3.count_ones() + (k as u32) < pl {
                continue;
            }
            self.slice(y0, y1, y3, &mut phi);
            let cx = self.ctx(y0, y1, y3);
            for y2 in 0..1u32 << k {
                let extra = if sym {
                    let ph = y3.count_ones() + y2.count_ones();
                    if ph < pl {
                        continue;
                    }
                    (ph > pl) as usize
                } else {
                    0
                };
                let f = phi[(!y2 as u64 & self.km) as usize] as u128;
                let (i, c) = self.c_fast(&cx, y2);
                hist[i * MAXC + c + extra] += f * f;
            }
        }
        self.combine(&hist)
    }

    /// Neighbours of the vertex set x in Q_d (not including x itself).
    #[inline]
    fn nbrs(&self, x: u64) -> u64 {
        let mut r = 0u64;
        for (j, &mj) in self.bitmask.iter().enumerate() {
            let s = 1u32 << j;
            r |= ((x & mj) << s) | ((x >> s) & mj);
        }
        let (k, n) = (self.k as u32, 4 * self.k as u32);
        let rot = if n == 64 {
            x.rotate_left(k) | x.rotate_right(k)
        } else {
            (x << k) | (x >> (n - k)) | (x >> k) | (x << (n - k))
        };
        (r | rot) & self.full
    }

    /// c(Y) = 3^isolated * 2^(components of size >= 2).
    #[inline]
    pub fn c(&self, y: u64) -> u128 {
        let iso = y & !self.nbrs(y);
        let mut rem = y & !iso;
        let mut comps = 0usize;
        while rem != 0 {
            let mut comp = rem & rem.wrapping_neg();
            loop {
                let nxt = (comp | self.nbrs(comp)) & y;
                if nxt == comp {
                    break;
                }
                comp = nxt;
            }
            rem &= !comp;
            comps += 1;
        }
        self.pow3[iso.count_ones() as usize] * self.pow2[comps]
    }

    /// The slice W2 -> f(W0, W1, W2, W3) for fixed W0, W1, W3 (all 2^k values of W2).
    pub fn slice(&self, y0: u64, y1: u64, y3: u64, phi: &mut [u64]) {
        let k = self.k;
        let (w0, w1, w3) = ((!y0 & self.km) as usize, (!y1 & self.km) as usize, (!y3 & self.km) as usize);
        let t0s = &self.lst[self.off[w0] as usize..self.off[w0 + 1] as usize];
        phi.iter_mut().for_each(|p| *p = 0);
        for &t2 in &self.im {
            let mut g = 0u64;
            for &t0 in t0s {
                let u = (t0 | t2) as usize;
                g += self.fm[w1 & !u] as u64 * self.fm[w3 & !u] as u64;
            }
            phi[t2 as usize] = g;
        }
        // subset-sum transform, branch-free so the inner loop vectorises
        for b in 0..k {
            let h = 1usize << b;
            for blk in phi.chunks_exact_mut(2 * h) {
                let (lo, hi) = blk.split_at_mut(h);
                for (x, y) in hi.iter_mut().zip(lo.iter()) {
                    *x += *y;
                }
            }
        }
    }

    /// sum over Y3, Y2 of c(Y) f(~Y)^2 for fixed Y_lo = Y0 | Y1 << k.
    /// `y3_range` restricts Y3 (for benchmarking); pass 0..2^k for the full sum.
    pub fn rep_sum(&self, y_lo: u64, y3_range: std::ops::Range<u64>) -> U256 {
        let k = self.k;
        let (y0, y1) = (y_lo & self.km, y_lo >> k & self.km);
        let mut phi = vec![0u64; 1usize << k];
        let mut acc = U256::default();
        for y3 in y3_range {
            self.slice(y0, y1, y3, &mut phi);
            let fixed = y0 | y1 << k | y3 << (3 * k);
            for y2 in 0..1u64 << k {
                let f = phi[(!y2 & self.km) as usize] as u128;
                let c = self.c(fixed | y2 << (2 * k));
                acc.add(U256::mul(c, f * f));
            }
        }
        acc
    }
}

/// d = 6 benchmark: random Y_lo, a few Y3 each; times the slice and the c(Y) loop separately and
/// spot-checks the slice against f_6 evaluated through Q_6 = Q_5 x K_2.
pub fn bench6(samples: usize, y3_per: u64) {
    let yg = YGroup::new(6);
    let fe5 = crate::FEval::new(5);
    let i5 = indep_sets(5);
    let mut rng = 0x2545F4914F6CDD1Du64;
    let mut next = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    // the 64-bit fibre layout is not exercised by the d <= 5 validation: check nbrs() and c()
    for v in 0..64u32 {
        let (f, x) = (v / 16, v % 16);
        let mut want = 0u64;
        for j in 0..4 {
            want |= 1 << (16 * f + (x ^ 1 << j));
        }
        want |= 1 << (16 * ((f + 1) % 4) + x) | 1 << (16 * ((f + 3) % 4) + x);
        assert_eq!(yg.nbrs(1 << v), want, "nbrs mismatch at vertex {v}");
    }
    // c(Y) against a brute-force count of pairs (S, U) of independent sets with S | U = Y,
    // on random small Y (the pairs are the 3^|Y| labellings that are proper)
    for trial in 0..300 {
        let mut y = 0u64;
        for _ in 0..(trial % 11 + 1) {
            y |= 1 << (next() % 64);
        }
        let vs: Vec<u32> = (0..64).filter(|&v| y >> v & 1 == 1).collect();
        let mut brute = 0u128;
        let mut lab = vec![0u8; vs.len()]; // 0 = S only, 1 = U only, 2 = both
        'outer: loop {
            let (mut s, mut u) = (0u64, 0u64);
            for (i, &v) in vs.iter().enumerate() {
                if lab[i] != 1 { s |= 1 << v; }
                if lab[i] != 0 { u |= 1 << v; }
            }
            if s & yg.nbrs(s) == 0 && u & yg.nbrs(u) == 0 { brute += 1; }
            for i in 0..vs.len() {
                lab[i] += 1;
                if lab[i] < 3 { continue 'outer; }
                lab[i] = 0;
            }
            break;
        }
        assert_eq!(brute, yg.c(y), "c mismatch for Y = {y:#x}");
    }
    println!("  d=6 layout: nbrs() correct for all 64 vertices; c(Y) = brute force on 300 random Y");
    let mut nfast = 0;
    for _ in 0..200 {
        let (y0, y1, y3) = (next() & 0xFFFF, next() & 0xFFFF, next() & 0xFFFF);
        let cx = yg.ctx(y0, y1, y3);
        for _ in 0..500 {
            let y2 = next() & 0xFFFF;
            let (i, c) = yg.c_fast(&cx, y2 as u32);
            assert_eq!(yg.pow3[i] * yg.pow2[c], yg.c(y0 | y1 << 16 | y2 << 32 | y3 << 48), "c_fast mismatch");
            nfast += 1;
        }
    }
    println!("  c_fast = c on {nfast} random d=6 terms");
    let (mut t_slice, mut t_c, mut nterms, mut checks) = (0f64, 0f64, 0u64, 0usize);
    let mut phi = vec![0u64; 1 << 16];
    let mut sink = U256::default();
    for _ in 0..samples {
        let y_lo = next() & 0xFFFF_FFFF;
        let (y0, y1) = (y_lo & 0xFFFF, y_lo >> 16);
        for _ in 0..y3_per {
            let y3 = next() & 0xFFFF;
            let t = Instant::now();
            yg.slice(y0, y1, y3, &mut phi);
            t_slice += t.elapsed().as_secs_f64();
            // spot check: lo Q_5 = fibres 0,1; hi Q_5 = fibre 3 (b = 0) and fibre 2 (b = 1)
            for _ in 0..4 {
                let y2 = next() & 0xFFFF;
                let wlo = (!(y0 | y1 << 16)) as u32;
                let whi = (!(y3 | y2 << 16)) as u32;
                let direct: u64 = i5.iter().filter(|&&a| a & !wlo == 0).map(|&a| fe5.f(whi & !a)).sum();
                assert_eq!(direct, phi[(!y2 & 0xFFFF) as usize], "slice mismatch");
                checks += 1;
            }
            let t = Instant::now();
            let cx = yg.ctx(y0, y1, y3);
            let mut hist = vec![0u128; MAXI * MAXC];
            for y2 in 0..1u32 << 16 {
                let f = phi[(!y2 & 0xFFFF) as usize] as u128;
                let (i, c) = yg.c_fast(&cx, y2);
                hist[i * MAXC + c] += f * f;
            }
            sink.add(yg.combine(&hist));
            t_c += t.elapsed().as_secs_f64();
            nterms += 1 << 16;
        }
    }
    let per_slice = t_slice / (samples as f64 * y3_per as f64);
    let per_term = t_c / nterms as f64;
    println!("bench d=6: {samples} random Y_lo x {y3_per} Y3, {checks} slice spot checks OK (sink {})", sink.lo % 10);
    println!("  slice (G + zeta over 2^16): {:.3} ms each = {:.1} ns per term", per_slice * 1e3, per_slice * 1e9 / 65536.0);
    println!("  c(Y) loop: {:.1} ns per term", per_term * 1e9);
    let terms = 1_228_158f64 * 4294967296.0;
    let per = per_term + per_slice / 65536.0;
    println!(
        "  full run: 1228158 Y_lo reps (A000616(5)) x 2^32 = {:.2e} terms -> {:.2e} core-seconds = {:.0} days on 12 threads",
        terms,
        terms * per,
        terms * per / 12.0 / 86400.0
    );
}

/// Full union-grouped sum for d <= 5 (orbit reps of Y_lo by union-find over all 2^(2k) masks).
pub fn run(d: usize) -> U256 {
    let t0 = Instant::now();
    let yg = YGroup::new(d);
    let nlo = 1u64 << (2 * yg.k);
    let all: Vec<u32> = (0..nlo as u32).collect();
    let reps = orbits(&all, &generators(d - 1, Group::Full));
    let (mut total, mut total_sym) = (U256::default(), U256::default());
    for &(r, w) in &reps {
        let fast = yg.rep_sum_fast(r as u64, 0..1u64 << yg.k, false);
        assert_eq!(fast, yg.rep_sum(r as u64, 0..1u64 << yg.k), "fast and slow rep sums differ");
        total.add(fast.scale(w));
        total_sym.add(yg.rep_sum_fast(r as u64, 0..1u64 << yg.k, true).scale(w));
    }
    assert_eq!(total, total_sym, "swap-symmetry weighting changes the total");
    eprintln!("ygroup d={d}: {} Y_lo orbit reps, plain and swap-weighted totals agree, {:.2}s", reps.len(), t0.elapsed().as_secs_f64());
    total
}

// ---------------------------------------------------------------------------------------------
// d = 6: orbit representatives of Y_lo and the CPU driver
// ---------------------------------------------------------------------------------------------

/// All 3840 elements of Aut(Q_5) as byte-table vertex permutations: v -> sigma(v) ^ a.
fn aut_q5() -> Vec<crate::VPerm> {
    fn perms(n: usize) -> Vec<Vec<usize>> {
        if n == 0 {
            return vec![vec![]];
        }
        let mut out = vec![];
        for p in perms(n - 1) {
            for i in 0..n {
                let mut q = p.clone();
                q.insert(i, n - 1);
                out.push(q);
            }
        }
        out
    }
    let mut g = vec![];
    for sigma in perms(5) {
        for a in 0..32usize {
            let p: Vec<usize> = (0..32usize)
                .map(|v| (0..5).fold(0, |w, i| w | (v >> i & 1) << sigma[i]) ^ a)
                .collect();
            g.push(crate::VPerm::new(&p));
        }
    }
    g
}

/// Orbit representatives (minimal elements) of all 2^32 subsets of V(Q_5) under Aut(Q_5), with
/// orbit sizes. Written as little-endian (u32 rep, u32 size) pairs.
pub fn enumerate_reps6(path: &std::path::Path) -> Vec<(u32, u32)> {
    let t0 = Instant::now();
    let g = aut_q5();
    assert_eq!(g.len(), 3840);
    let mut seen = vec![0u64; 1 << 26]; // one bit per mask, 512 MB
    let mut reps = vec![];
    let mut total = 0u64;
    for word in 0..(1usize << 26) {
        while seen[word] != u64::MAX {
            let bit = (!seen[word]).trailing_zeros() as usize;
            let x = (word << 6 | bit) as u32;
            let mut size = 0u32;
            for e in &g {
                let y = e.apply(x) as usize;
                let (w, b) = (y >> 6, y & 63);
                if seen[w] >> b & 1 == 0 {
                    seen[w] |= 1 << b;
                    size += 1;
                }
            }
            reps.push((x, size));
            total += size as u64;
        }
        if word & 0xFFFFF == 0xFFFFF {
            eprintln!("  enumerate: {:.0}% ({} reps, {:.0}s)", 100.0 * (word + 1) as f64 / (1u64 << 26) as f64, reps.len(), t0.elapsed().as_secs_f64());
        }
    }
    assert_eq!(total, 1u64 << 32, "orbit sizes must cover every subset");
    let mut bytes = Vec::with_capacity(reps.len() * 8);
    for &(x, w) in &reps {
        bytes.extend_from_slice(&x.to_le_bytes());
        bytes.extend_from_slice(&w.to_le_bytes());
    }
    std::fs::write(path, bytes).unwrap();
    eprintln!("enumerate: {} orbit reps (A000616(5) = 1228158), {:.0}s -> {}", reps.len(), t0.elapsed().as_secs_f64(), path.display());
    reps
}

pub fn load_reps6(path: &std::path::Path) -> Vec<(u32, u32)> {
    let b = std::fs::read(path).expect("run `hypercube-indep reps6` first");
    b.chunks_exact(8)
        .map(|c| (u32::from_le_bytes(c[0..4].try_into().unwrap()), u32::from_le_bytes(c[4..8].try_into().unwrap())))
        .collect()
}

/// Unweighted rep sums for chosen rep indices (reference values for checking the GPU),
/// one rep per thread.
pub fn cpu_reps6(reps: &[(u32, u32)], idx: &[usize], sym: bool) -> Vec<U256> {
    let yg = YGroup::new(6);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let out = std::sync::Mutex::new(vec![U256::default(); idx.len()]);
    std::thread::scope(|sc| {
        for _ in 0..crate::num_threads() {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= idx.len() {
                    break;
                }
                let t = Instant::now();
                let r = yg.rep_sum_fast(reps[idx[i]].0 as u64, 0..1 << 16, sym);
                eprintln!("  rep #{} (Y_lo = {:#010x}, orbit {}): {}  [{:.1}s]", idx[i], reps[idx[i]].0, reps[idx[i]].1, r.to_decimal(), t.elapsed().as_secs_f64());
                out.lock().unwrap()[i] = r;
            });
        }
    });
    out.into_inner().unwrap()
}
