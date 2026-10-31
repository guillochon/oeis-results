//! A000080 (minimal triangle graphs) by species theory, modulo a prime.
//!
//! See ../README.md for the derivation.  In short: a minimal connected triple system is a tree of points and
//! units, a unit being a lone triple E3(X) or a 2-connected loopless multigraph block B(X,Y) whose edges (Y)
//! each carry one pendant "arm" point.  B is extracted from the connected edge-labelled multigraph species C
//! through the block decomposition, all evaluated by plethystic substitution into one-variable series:
//!
//!   G(S,W)  = sum_lambda p_lambda(S)/z_lambda * exp(sum_d o_d(lambda) L(x^d)),  L = sum_j W(x^j)/j
//!   C^.(S,W) = R   (defines S),   B*_X(R,R) = PLog(R/S),   B*_Y(R,R) = (dG/dq1)/G
//!   R = x PExp(E2(R) + B_X(R,R) + B_Y(R,R)),   H = R + U(R) - R U'(R)
//!
//! Usage: mtg N PRIME_INDEX THREADS   -> prints "n a(n) mod p" for n = 3..N.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

const PRIMES: [u64; 5] = [
    4611686018427387847,
    4611686018427387817,
    4611686018427387787,
    4611686018427387761,
    4611686018427387751,
];

#[derive(Clone, Copy)]
struct Fp {
    p: u64,
}

impl Fp {
    #[inline]
    fn add(self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }
    #[inline]
    fn sub(self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }
    #[inline]
    fn mul(self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn pow(self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    fn inv(self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
    fn neg(self, a: u64) -> u64 {
        if a == 0 { 0 } else { self.p - a }
    }
    #[inline]
    fn reduce(self, acc: u128) -> u64 {
        (acc % self.p as u128) as u64
    }
}

/// Dot-product accumulator: products are < 2^124, so reduce before the sum can pass 2^127.
macro_rules! acc_add {
    ($f:expr, $acc:ident, $x:expr) => {
        $acc += $x;
        if $acc >> 126 != 0 {
            $acc %= $f.p as u128;
        }
    };
}

struct Ctx {
    f: Fp,
    inv: Vec<u64>,
    mu: Vec<i64>,
}

impl Ctx {
    fn new(p: u64, n: usize) -> Ctx {
        let f = Fp { p };
        let mut inv = vec![0u64; n + 2];
        for i in 1..n + 2 {
            inv[i] = f.inv(i as u64);
        }
        let mut mu = vec![0i64; n + 2];
        for k in 1..n + 2 {
            mu[k] = mobius(k);
        }
        Ctx { f, inv, mu }
    }

    fn mul(&self, a: &[u64], b: &[u64], len: usize) -> Vec<u64> {
        let f = self.f;
        let mut c = vec![0u64; len];
        for n in 0..len {
            let mut acc: u128 = 0;
            for i in 0..=n.min(a.len().saturating_sub(1)) {
                if n - i < b.len() {
                    acc_add!(f, acc, a[i] as u128 * b[n - i] as u128);
                }
            }
            c[n] = f.reduce(acc);
        }
        c
    }

    fn exp(&self, a: &[u64], len: usize) -> Vec<u64> {
        let f = self.f;
        let ka: Vec<u64> = (0..len).map(|k| if k < a.len() { f.mul(k as u64 % f.p, a[k]) } else { 0 }).collect();
        let mut e = vec![0u64; len];
        e[0] = 1;
        for n in 1..len {
            let mut acc: u128 = 0;
            for k in 1..=n {
                acc_add!(f, acc, ka[k] as u128 * e[n - k] as u128);
            }
            e[n] = f.mul(f.reduce(acc), self.inv[n]);
        }
        e
    }

    fn log(&self, a: &[u64], len: usize) -> Vec<u64> {
        // a[0] == 1
        let f = self.f;
        let mut kl = vec![0u64; len]; // k * l[k]
        let mut l = vec![0u64; len];
        for n in 1..len {
            let mut acc: u128 = 0;
            for k in 1..n {
                acc_add!(f, acc, kl[k] as u128 * a[n - k] as u128);
            }
            let v = f.sub(f.mul(n as u64, a[n]), f.reduce(acc));
            kl[n] = v;
            l[n] = f.mul(v, self.inv[n]);
        }
        l
    }

    fn inverse(&self, a: &[u64], len: usize) -> Vec<u64> {
        let f = self.f;
        let i0 = f.inv(a[0]);
        let mut b = vec![0u64; len];
        b[0] = i0;
        for n in 1..len {
            let mut acc: u128 = 0;
            for k in 1..=n.min(a.len() - 1) {
                acc_add!(f, acc, a[k] as u128 * b[n - k] as u128);
            }
            b[n] = f.mul(f.neg(f.reduce(acc)), i0);
        }
        b
    }

    /// a(x^k) truncated to len
    fn at(&self, a: &[u64], k: usize, len: usize) -> Vec<u64> {
        let mut c = vec![0u64; len];
        let mut i = 0;
        while i * k < len && i < a.len() {
            c[i * k] = a[i];
            i += 1;
        }
        c
    }

    fn lincomb_at(&self, out: &mut [u64], a: &[u64], k: usize, coef: u64) {
        let f = self.f;
        let mut i = 1;
        while i * k < out.len() && i < a.len() {
            out[i * k] = f.add(out[i * k], f.mul(coef, a[i]));
            i += 1;
        }
    }

    fn mu_over_k(&self, k: usize) -> u64 {
        let f = self.f;
        match self.mu[k] {
            0 => 0,
            1 => self.inv[k],
            _ => f.neg(self.inv[k]),
        }
    }

    /// PExp(a) = exp(sum_k a(x^k)/k), a[0] == 0
    fn pexp(&self, a: &[u64], len: usize) -> Vec<u64> {
        let mut s = vec![0u64; len];
        for k in 1..len {
            self.lincomb_at(&mut s, a, k, self.inv[k]);
        }
        self.exp(&s, len)
    }

    /// PLog(a) = sum_k mu(k)/k log(a)(x^k), a[0] == 1
    fn plog(&self, a: &[u64], len: usize) -> Vec<u64> {
        let l = self.log(a, len);
        let mut s = vec![0u64; len];
        for k in 1..len {
            let c = self.mu_over_k(k);
            if c != 0 {
                self.lincomb_at(&mut s, &l, k, c);
            }
        }
        s
    }

    fn e2(&self, r: &[u64], len: usize) -> Vec<u64> {
        let f = self.f;
        let rr = self.mul(r, r, len);
        let r2 = self.at(r, 2, len);
        (0..len).map(|i| f.mul(f.add(rr[i], r2[i]), self.inv[2])).collect()
    }

    fn e3(&self, r: &[u64], len: usize) -> Vec<u64> {
        let f = self.f;
        let rr = self.mul(r, r, len);
        let rrr = self.mul(&rr, r, len);
        let r2 = self.at(r, 2, len);
        let rr2 = self.mul(r, &r2, len);
        let r3 = self.at(r, 3, len);
        (0..len)
            .map(|i| {
                let v = f.add(f.add(rrr[i], f.mul(3, rr2[i])), f.mul(2, r3[i]));
                f.mul(v, self.inv[6])
            })
            .collect()
    }
}

fn mobius(mut n: usize) -> i64 {
    let mut r = 1;
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return 0;
            }
            r = -r;
        }
        p += 1;
    }
    if n > 1 { -r } else { r }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Per-thread DFS over partitions lambda of n <= t.
struct Walker<'a> {
    c: &'a Ctx,
    t: usize,
    s: &'a [u64],
    l: &'a [u64],
    parts: Vec<(usize, usize)>, // (part, multiplicity), parts decreasing
    prod: Vec<Vec<u64>>,        // per depth
    ex: Vec<Vec<u64>>,          // per depth
    g: Vec<u64>,
    gp: Vec<u64>,
    gq: Vec<u64>,
}

impl<'a> Walker<'a> {
    fn new(c: &'a Ctx, t: usize, s: &'a [u64], l: &'a [u64]) -> Self {
        let mut prod = vec![vec![0u64; t + 1]; t + 2];
        prod[0][0] = 1;
        Walker {
            c,
            t,
            s,
            l,
            parts: Vec::new(),
            prod,
            ex: vec![vec![0u64; t + 1]; t + 2],
            g: vec![0; t + 1],
            gp: vec![0; t + 1],
            gq: vec![0; t + 1],
        }
    }

    /// Push part a onto the partition at depth d (total n before), filling prod[d+1], ex[d+1].
    /// Returns (new n, new zinv).
    fn push(&mut self, d: usize, n: usize, zinv: u64, a: usize) -> (usize, u64) {
        let f = self.c.f;
        let t = self.t;
        let nn = n + a;
        let lim = t - nn; // degrees of ex / offsets of prod that still matter
        // multiplicity
        let m = match self.parts.last_mut() {
            Some((b, m)) if *b == a => {
                *m += 1;
                *m
            }
            _ => {
                self.parts.push((a, 1));
                1
            }
        };
        let zinv2 = f.mul(f.mul(zinv, self.c.inv[a]), self.c.inv[m]);
        // prod[d+1] = prod[d] * S(x^a), support [nn, t]
        {
            let (lo, hi) = self.prod.split_at_mut(d + 1);
            let src = &lo[d];
            let dst = &mut hi[0];
            for o in 0..=lim {
                // dst[nn+o] = sum_{j>=1, a*j <= a + o} src[nn + o - a*j] * s[j]
                let mut acc: u128 = 0;
                let mut j = 1;
                while a * j <= a + o {
                    let idx = nn + o - a * j;
                    acc_add!(f, acc, src[idx] as u128 * self.s[j] as u128);
                    j += 1;
                }
                dst[nn + o] = f.reduce(acc);
            }
        }
        // ex[d+1] = ex[d] + pair orbits created by the new cycle
        {
            let (lo, hi) = self.ex.split_at_mut(d + 1);
            let dst = &mut hi[0];
            dst[..=lim].copy_from_slice(&lo[d][..=lim]);
            let mut addorb = |len: usize, cnt: usize| {
                if cnt == 0 || len > lim {
                    return;
                }
                let cm = cnt as u64 % f.p;
                let mut i = 1;
                while i * len <= lim {
                    dst[i * len] = f.add(dst[i * len], f.mul(cm, self.l[i]));
                    i += 1;
                }
            };
            if a % 2 == 1 {
                addorb(a, (a - 1) / 2);
            } else {
                addorb(a, a / 2 - 1);
                addorb(a / 2, 1);
            }
            for &(b, mb) in self.parts.iter() {
                let k = if b == a { mb - 1 } else { mb };
                if k > 0 {
                    let g = gcd(a, b);
                    addorb(a / g * b, k * g);
                }
            }
        }
        (nn, zinv2)
    }

    fn pop(&mut self) {
        let last = self.parts.last_mut().unwrap();
        if last.1 == 1 {
            self.parts.pop();
        } else {
            last.1 -= 1;
        }
    }

    /// Add the term of the current partition (at depth d, total n).
    fn term(&mut self, d: usize, n: usize, zinv: u64) {
        let f = self.c.f;
        let lim = self.t - n;
        let e = self.c.exp(&self.ex[d][..=lim], lim + 1);
        let pr = &self.prod[d];
        let m1 = self.parts.iter().find(|p| p.0 == 1).map_or(0, |p| p.1);
        let m2 = self.parts.iter().find(|p| p.0 == 2).map_or(0, |p| p.1);
        let w1 = m1 as u64;
        let wq = (m1 * m1.saturating_sub(1) / 2 + m2) as u64;
        for s in 0..=lim {
            let mut acc: u128 = 0;
            for j in 0..=s {
                acc_add!(f, acc, pr[n + j] as u128 * e[s - j] as u128);
            }
            let v = f.mul(f.reduce(acc), zinv);
            let i = n + s;
            self.g[i] = f.add(self.g[i], v);
            if w1 != 0 {
                self.gp[i] = f.add(self.gp[i], f.mul(v, w1));
            }
            if wq != 0 {
                self.gq[i] = f.add(self.gq[i], f.mul(v, wq));
            }
        }
    }

    fn rec(&mut self, d: usize, n: usize, zinv: u64, maxpart: usize) {
        for a in 1..=maxpart.min(self.t - n) {
            let (nn, z2) = self.push(d, n, zinv, a);
            self.term(d + 1, nn, z2);
            self.rec(d + 1, nn, z2, a);
            self.pop();
        }
    }
}

/// G, p1 dG/dp1, dG/dq1 at p_k = S(x^k), q_k = W(x^k), truncated at degree t.
fn graph_series(c: &Ctx, s: &[u64], w: &[u64], t: usize, threads: usize) -> (Vec<u64>, Vec<u64>, Vec<u64>) {
    let f = c.f;
    let len = t + 1;
    let mut l = vec![0u64; len];
    for j in 1..len {
        c.lincomb_at(&mut l, &w[..len.min(w.len())], j, c.inv[j]);
    }
    let s = &s[..len];
    // tasks: (a1, a2) with a2 = 0 meaning the single-part node
    let mut tasks = Vec::new();
    for a1 in 1..=t {
        tasks.push((a1, 0));
        for a2 in 1..=a1.min(t - a1) {
            tasks.push((a1, a2));
        }
    }
    // heaviest first: small parts leave the most room
    tasks.sort_by_key(|&(a1, a2)| a1 + a2);
    let next = AtomicUsize::new(0);
    let out = Mutex::new((vec![0u64; len], vec![0u64; len], vec![0u64; len]));
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| {
                let mut wk = Walker::new(c, t, s, &l);
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= tasks.len() {
                        break;
                    }
                    let (a1, a2) = tasks[i];
                    let (n1, z1) = wk.push(0, 0, 1, a1);
                    if a2 == 0 {
                        wk.term(1, n1, z1);
                    } else {
                        let (n2, z2) = wk.push(1, n1, z1, a2);
                        wk.term(2, n2, z2);
                        wk.rec(2, n2, z2, a2);
                        wk.pop();
                    }
                    wk.pop();
                }
                let mut o = out.lock().unwrap();
                for i in 0..len {
                    o.0[i] = f.add(o.0[i], wk.g[i]);
                    o.1[i] = f.add(o.1[i], wk.gp[i]);
                    o.2[i] = f.add(o.2[i], wk.gq[i]);
                }
            });
        }
    });
    let (mut g, gp, gq) = out.into_inner().unwrap();
    g[0] = f.add(g[0], 1); // empty graph
    (g, gp, gq)
}

/// U'(R) = E2(R) + B_X(R,R) + B_Y(R,R) and the pieces needed for U, from G evaluated at (S,R).
struct Pieces {
    bxs: Vec<u64>,
    cy: Vec<u64>,
}

fn pieces(c: &Ctx, r: &[u64], s: &[u64], g: &[u64], gq: &[u64], len: usize) -> Pieces {
    // B*_X(R,R) = PLog((R/x)/(S/x)), needs R, S to degree len
    let q: Vec<u64> = (0..len).map(|i| r[i + 1]).collect();
    let sx: Vec<u64> = (0..len).map(|i| s[i + 1]).collect();
    let ratio = c.mul(&q, &c.inverse(&sx, len), len);
    let bxs = c.plog(&ratio, len);
    let gi = c.inverse(g, len);
    let cy = c.mul(gq, &gi, len);
    Pieces { bxs, cy }
}

fn uprime(c: &Ctx, r: &[u64], p: &Pieces, len: usize) -> Vec<u64> {
    let f = c.f;
    let rr = c.mul(r, r, len);
    let e2 = c.e2(r, len);
    // E2 + (B*_X - R^2) + (B*_Y - E2) = B*_X - R^2 + B*_Y
    let _ = e2;
    (0..len).map(|i| f.add(f.sub(p.bxs[i], rr[i]), p.cy[i])).collect()
}

/// Windows 11 keeps background console jobs on the E-cores (EcoQoS); opt out of that throttling, and run at
/// below-normal priority so the machine stays responsive.
#[cfg(windows)]
fn be_a_good_citizen() {
    #[repr(C)]
    struct ProcessPowerThrottlingState {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
        fn SetPriorityClass(h: isize, class: u32) -> i32;
    }
    const PROCESS_POWER_THROTTLING: i32 = 4;
    const EXECUTION_SPEED: u32 = 0x1;
    const BELOW_NORMAL_PRIORITY_CLASS: u32 = 0x4000;
    let st = ProcessPowerThrottlingState { version: 1, control_mask: EXECUTION_SPEED, state_mask: 0 };
    unsafe {
        if SetProcessInformation(
            GetCurrentProcess(),
            PROCESS_POWER_THROTTLING,
            &st as *const _ as *const std::ffi::c_void,
            std::mem::size_of::<ProcessPowerThrottlingState>() as u32,
        ) == 0
        {
            eprintln!("note: could not disable Windows power throttling (EcoQoS)");
        }
        SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
    }
}

#[cfg(not(windows))]
fn be_a_good_citizen() {}

fn main() {
    be_a_good_citizen();
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).map_or(20, |s| s.parse().unwrap());
    let pi: usize = args.get(2).map_or(0, |s| s.parse().unwrap());
    let threads: usize = args.get(3).map_or(6, |s| s.parse().unwrap());
    let p = PRIMES[pi];
    let c = Ctx::new(p, n + 2);
    let f = c.f;
    let t0 = Instant::now();
    // R, S with one spare slot (index n+1) used by the ratio R/S
    let mut r = vec![0u64; n + 2];
    let mut s = vec![0u64; n + 2];
    r[1] = 1;
    s[1] = 1;
    for k in 2..=n {
        let (g, gp, gq) = graph_series(&c, &s, &r, k, threads);
        let len = k + 1;
        let cdot = c.mul(&gp, &c.inverse(&g, len), len);
        let dk = cdot[k]; // = S_k + D_k with S_k = 0
        // ratio at order k-1 needs R_k - S_k = dk
        let mut rt = r.clone();
        rt[k] = dk;
        let pc = pieces(&c, &rt, &s, &g[..k], &gq[..k], k);
        let ud = uprime(&c, &r, &pc, k);
        let pe = c.pexp(&ud, k);
        r[k] = pe[k - 1];
        s[k] = f.sub(r[k], dk);
        eprintln!("k={k} {:.1}s", t0.elapsed().as_secs_f64());
    }
    // final assembly at degree n
    let len = n + 1;
    let (g, _gp, gq) = graph_series(&c, &s, &r, n, threads);
    // need R_{n+1} - S_{n+1} only at ratio degree n, which only enters B*_X at degree n (not used in H)
    let pc = pieces(&c, &r, &s, &g, &gq, len);
    let lg = c.log(&g, len);
    let mut cc = vec![0u64; len];
    for k in 1..len {
        let m = c.mu_over_k(k);
        if m != 0 {
            c.lincomb_at(&mut cc, &lg, k, m);
        }
    }
    let rr = &r[..len];
    let rb = c.mul(rr, &pc.bxs, len);
    let e2 = c.e2(rr, len);
    let e2r = c.mul(&e2, rr, len);
    // B(R,R) = C - R + R B*_X - E2(R) R ; U = E3(R) + B
    let e3 = c.e3(rr, len);
    let u: Vec<u64> = (0..len)
        .map(|i| f.add(e3[i], f.sub(f.add(f.sub(cc[i], rr[i]), rb[i]), e2r[i])))
        .collect();
    let ud = uprime(&c, rr, &pc, len);
    let rud = c.mul(rr, &ud, len);
    for i in 3..len {
        let h = f.sub(f.add(rr[i], u[i]), rud[i]);
        println!("{i} {h}");
    }
    eprintln!("done {:.1}s", t0.elapsed().as_secs_f64());
}
