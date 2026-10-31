//! A000080 and A048781 together: the main.rs formula with every triple weighted by y, evaluated over dual
//! numbers y = 1 + eps (eps^2 = 0).  The constant part counts minimal triangle graphs (A000080), the eps part
//! counts their triples (A048781).  Plethystic substitution sends y -> y^k = 1 + k eps.
//!
//!   W = yR,  U' = B*_X - y R^2 + y B*_Y,  U = y E3(R) + B*(R,W) - y E2(R) R,  H = R + U - R U'
//!
//! Usage: mtg2 N PRIME_INDEX THREADS   -> prints "n A000080(n) A048781(n)" mod p for n = 3..N.

use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

const PRIMES: [u64; 5] = [
    4611686018427387847,
    4611686018427387817,
    4611686018427387787,
    4611686018427387761,
    4611686018427387751,
];

#[derive(Clone, Copy, Default, PartialEq)]
struct D {
    a: u64,
    b: u64,
}

const ZERO: D = D { a: 0, b: 0 };
const ONE: D = D { a: 1, b: 0 };

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
    // dual-number operations
    #[inline]
    fn dadd(self, x: D, y: D) -> D {
        D { a: self.add(x.a, y.a), b: self.add(x.b, y.b) }
    }
    #[inline]
    fn dsub(self, x: D, y: D) -> D {
        D { a: self.sub(x.a, y.a), b: self.sub(x.b, y.b) }
    }
    #[inline]
    fn dmul(self, x: D, y: D) -> D {
        D { a: self.mul(x.a, y.a), b: self.add(self.mul(x.a, y.b), self.mul(x.b, y.a)) }
    }
    #[inline]
    fn dscale(self, x: D, c: u64) -> D {
        D { a: self.mul(x.a, c), b: self.mul(x.b, c) }
    }
    /// x with y -> y^k, i.e. eps -> k eps
    #[inline]
    fn dpow_y(self, x: D, k: u64) -> D {
        D { a: x.a, b: self.mul(x.b, k % self.p) }
    }
    fn dinv(self, x: D) -> D {
        let ia = self.inv(x.a);
        D { a: ia, b: self.neg(self.mul(x.b, self.mul(ia, ia))) }
    }
}

/// Dot products of duals: accumulate the a and b parts in u128, reducing before 2^127.
struct Acc {
    a: u128,
    b: u128,
}

impl Acc {
    #[inline]
    fn new() -> Acc {
        Acc { a: 0, b: 0 }
    }
    #[inline]
    fn add(&mut self, p: u64, x: D, y: D) {
        self.a += x.a as u128 * y.a as u128;
        self.b += x.a as u128 * y.b as u128;
        self.b += x.b as u128 * y.a as u128;
        if (self.a | self.b) >> 126 != 0 {
            self.a %= p as u128;
            self.b %= p as u128;
        }
    }
    #[inline]
    fn get(&self, f: Fp) -> D {
        D { a: f.reduce(self.a), b: f.reduce(self.b) }
    }
}

type Ser = Vec<D>;

struct Ctx {
    f: Fp,
    inv: Vec<u64>,
    mu: Vec<i64>,
}

impl Ctx {
    fn new(p: u64, n: usize) -> Ctx {
        let f = Fp { p };
        let inv = (0..n + 2).map(|i| if i == 0 { 0 } else { f.inv(i as u64) }).collect();
        let mu = (0..n + 2).map(mobius).collect();
        Ctx { f, inv, mu }
    }

    fn mul(&self, a: &[D], b: &[D], len: usize) -> Ser {
        let p = self.f.p;
        (0..len)
            .map(|n| {
                let mut acc = Acc::new();
                for i in 0..=n.min(a.len().saturating_sub(1)) {
                    if n - i < b.len() {
                        acc.add(p, a[i], b[n - i]);
                    }
                }
                acc.get(self.f)
            })
            .collect()
    }

    fn exp(&self, a: &[D], len: usize) -> Ser {
        let f = self.f;
        let ka: Ser = (0..len).map(|k| if k < a.len() { f.dscale(a[k], k as u64) } else { ZERO }).collect();
        let mut e = vec![ZERO; len];
        e[0] = ONE;
        for n in 1..len {
            let mut acc = Acc::new();
            for k in 1..=n {
                acc.add(f.p, ka[k], e[n - k]);
            }
            e[n] = f.dscale(acc.get(f), self.inv[n]);
        }
        e
    }

    fn log(&self, a: &[D], len: usize) -> Ser {
        // a[0] == 1
        let f = self.f;
        let mut kl = vec![ZERO; len];
        let mut l = vec![ZERO; len];
        for n in 1..len {
            let mut acc = Acc::new();
            for k in 1..n {
                acc.add(f.p, kl[k], a[n - k]);
            }
            let v = f.dsub(f.dscale(a[n], n as u64), acc.get(f));
            kl[n] = v;
            l[n] = f.dscale(v, self.inv[n]);
        }
        l
    }

    fn inverse(&self, a: &[D], len: usize) -> Ser {
        let f = self.f;
        let i0 = f.dinv(a[0]);
        let mut b = vec![ZERO; len];
        b[0] = i0;
        for n in 1..len {
            let mut acc = Acc::new();
            for k in 1..=n.min(a.len() - 1) {
                acc.add(f.p, a[k], b[n - k]);
            }
            b[n] = f.dmul(f.dsub(ZERO, acc.get(f)), i0);
        }
        b
    }

    /// a(x^k, y^k) truncated to len
    fn at(&self, a: &[D], k: usize, len: usize) -> Ser {
        let mut c = vec![ZERO; len];
        let mut i = 0;
        while i * k < len && i < a.len() {
            c[i * k] = self.f.dpow_y(a[i], k as u64);
            i += 1;
        }
        c
    }

    /// out += coef * a(x^k, y^k), skipping the constant term
    fn lincomb_at(&self, out: &mut [D], a: &[D], k: usize, coef: u64) {
        let f = self.f;
        let mut i = 1;
        while i * k < out.len() && i < a.len() {
            out[i * k] = f.dadd(out[i * k], f.dscale(f.dpow_y(a[i], k as u64), coef));
            i += 1;
        }
    }

    fn mu_over_k(&self, k: usize) -> u64 {
        match self.mu[k] {
            0 => 0,
            1 => self.inv[k],
            _ => self.f.neg(self.inv[k]),
        }
    }

    fn pexp(&self, a: &[D], len: usize) -> Ser {
        let mut s = vec![ZERO; len];
        for k in 1..len {
            self.lincomb_at(&mut s, a, k, self.inv[k]);
        }
        self.exp(&s, len)
    }

    fn plog(&self, a: &[D], len: usize) -> Ser {
        let l = self.log(a, len);
        let mut s = vec![ZERO; len];
        for k in 1..len {
            let c = self.mu_over_k(k);
            if c != 0 {
                self.lincomb_at(&mut s, &l, k, c);
            }
        }
        s
    }

    fn e2(&self, r: &[D], len: usize) -> Ser {
        let f = self.f;
        let rr = self.mul(r, r, len);
        let r2 = self.at(r, 2, len);
        (0..len).map(|i| f.dscale(f.dadd(rr[i], r2[i]), self.inv[2])).collect()
    }

    fn e3(&self, r: &[D], len: usize) -> Ser {
        let f = self.f;
        let rr = self.mul(r, r, len);
        let rrr = self.mul(&rr, r, len);
        let r2 = self.at(r, 2, len);
        let rr2 = self.mul(r, &r2, len);
        let r3 = self.at(r, 3, len);
        (0..len)
            .map(|i| f.dscale(f.dadd(f.dadd(rrr[i], f.dscale(rr2[i], 3)), f.dscale(r3[i], 2)), self.inv[6]))
            .collect()
    }

    fn times_y(&self, a: &[D]) -> Ser {
        // y = 1 + eps
        a.iter().map(|&x| D { a: x.a, b: self.f.add(x.a, x.b) }).collect()
    }
}

fn mobius(mut n: usize) -> i64 {
    if n == 0 {
        return 0;
    }
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

struct Walker<'a> {
    c: &'a Ctx,
    t: usize,
    s: &'a [D],
    l: &'a [D],
    parts: Vec<(usize, usize)>,
    prod: Vec<Ser>,
    ex: Vec<Ser>,
    g: Ser,
    gp: Ser,
    gq: Ser,
}

impl<'a> Walker<'a> {
    fn new(c: &'a Ctx, t: usize, s: &'a [D], l: &'a [D]) -> Self {
        let mut prod = vec![vec![ZERO; t + 1]; t + 2];
        prod[0][0] = ONE;
        Walker {
            c,
            t,
            s,
            l,
            parts: Vec::new(),
            prod,
            ex: vec![vec![ZERO; t + 1]; t + 2],
            g: vec![ZERO; t + 1],
            gp: vec![ZERO; t + 1],
            gq: vec![ZERO; t + 1],
        }
    }

    fn push(&mut self, d: usize, n: usize, zinv: u64, a: usize) -> (usize, u64) {
        let f = self.c.f;
        let t = self.t;
        let nn = n + a;
        let lim = t - nn;
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
        {
            let (lo, hi) = self.prod.split_at_mut(d + 1);
            let src = &lo[d];
            let dst = &mut hi[0];
            for o in 0..=lim {
                let mut acc = Acc::new();
                let mut j = 1;
                while a * j <= a + o {
                    acc.add(f.p, src[nn + o - a * j], f.dpow_y(self.s[j], a as u64));
                    j += 1;
                }
                dst[nn + o] = acc.get(f);
            }
        }
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
                    dst[i * len] = f.dadd(dst[i * len], f.dscale(f.dpow_y(self.l[i], len as u64), cm));
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
            let mut acc = Acc::new();
            for j in 0..=s {
                acc.add(f.p, pr[n + j], e[s - j]);
            }
            let v = f.dscale(acc.get(f), zinv);
            let i = n + s;
            self.g[i] = f.dadd(self.g[i], v);
            if w1 != 0 {
                self.gp[i] = f.dadd(self.gp[i], f.dscale(v, w1));
            }
            if wq != 0 {
                self.gq[i] = f.dadd(self.gq[i], f.dscale(v, wq));
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

fn graph_series(c: &Ctx, s: &[D], w: &[D], t: usize, threads: usize) -> (Ser, Ser, Ser) {
    let f = c.f;
    let len = t + 1;
    let mut l = vec![ZERO; len];
    for j in 1..len {
        c.lincomb_at(&mut l, &w[..len.min(w.len())], j, c.inv[j]);
    }
    let s = &s[..len];
    let mut tasks = Vec::new();
    for a1 in 1..=t {
        tasks.push((a1, 0));
        for a2 in 1..=a1.min(t - a1) {
            tasks.push((a1, a2));
        }
    }
    tasks.sort_by_key(|&(a1, a2)| a1 + a2);
    let next = AtomicUsize::new(0);
    let out = Mutex::new((vec![ZERO; len], vec![ZERO; len], vec![ZERO; len]));
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
                    o.0[i] = f.dadd(o.0[i], wk.g[i]);
                    o.1[i] = f.dadd(o.1[i], wk.gp[i]);
                    o.2[i] = f.dadd(o.2[i], wk.gq[i]);
                }
            });
        }
    });
    let (mut g, gp, gq) = out.into_inner().unwrap();
    g[0] = f.dadd(g[0], ONE);
    (g, gp, gq)
}

/// B*_X(R,W) = PLog(R/S) and B*_Y(R,W) = (dG/dq1)/G, to degree len-1
fn pieces(c: &Ctx, r: &[D], s: &[D], g: &[D], gq: &[D], len: usize) -> (Ser, Ser) {
    let q: Ser = (0..len).map(|i| r[i + 1]).collect();
    let sx: Ser = (0..len).map(|i| s[i + 1]).collect();
    let ratio = c.mul(&q, &c.inverse(&sx, len), len);
    let bxs = c.plog(&ratio, len);
    let cy = c.mul(gq, &c.inverse(g, len), len);
    (bxs, cy)
}

/// U' = B*_X - y R^2 + y B*_Y
fn uprime(c: &Ctx, r: &[D], bxs: &[D], cy: &[D], len: usize) -> Ser {
    let f = c.f;
    let rr = c.times_y(&c.mul(r, r, len));
    let ycy = c.times_y(&cy[..len]);
    (0..len).map(|i| f.dadd(f.dsub(bxs[i], rr[i]), ycy[i])).collect()
}

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
    let st = ProcessPowerThrottlingState { version: 1, control_mask: 1, state_mask: 0 };
    unsafe {
        if SetProcessInformation(
            GetCurrentProcess(),
            4,
            &st as *const _ as *const std::ffi::c_void,
            std::mem::size_of::<ProcessPowerThrottlingState>() as u32,
        ) == 0
        {
            eprintln!("note: could not disable Windows power throttling (EcoQoS)");
        }
        SetPriorityClass(GetCurrentProcess(), 0x4000);
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
    let c = Ctx::new(PRIMES[pi], n + 2);
    let f = c.f;
    let t0 = Instant::now();
    let mut r = vec![ZERO; n + 2];
    let mut s = vec![ZERO; n + 2];
    r[1] = ONE;
    s[1] = ONE;
    for k in 2..=n {
        let w = c.times_y(&r);
        let (g, gp, gq) = graph_series(&c, &s, &w, k, threads);
        let len = k + 1;
        let cdot = c.mul(&gp, &c.inverse(&g, len), len);
        let dk = cdot[k];
        let mut rt = r.clone();
        rt[k] = dk;
        let (bxs, cy) = pieces(&c, &rt, &s, &g[..k], &gq[..k], k);
        let ud = uprime(&c, &r, &bxs, &cy, k);
        let pe = c.pexp(&ud, k);
        r[k] = pe[k - 1];
        s[k] = f.dsub(r[k], dk);
        eprintln!("k={k} {:.1}s", t0.elapsed().as_secs_f64());
    }
    let len = n + 1;
    let w = c.times_y(&r);
    let (g, _gp, gq) = graph_series(&c, &s, &w, n, threads);
    let (bxs, cy) = pieces(&c, &r, &s, &g, &gq, len);
    let lg = c.log(&g, len);
    let mut cc = vec![ZERO; len];
    for k in 1..len {
        let m = c.mu_over_k(k);
        if m != 0 {
            c.lincomb_at(&mut cc, &lg, k, m);
        }
    }
    let rr = &r[..len];
    let rb = c.mul(rr, &bxs, len);
    let e2r = c.times_y(&c.mul(&c.e2(rr, len), rr, len));
    let ye3 = c.times_y(&c.e3(rr, len));
    // U = y E3(R) + (C - R + R B*_X) - y E2(R) R
    let u: Ser = (0..len).map(|i| f.dadd(ye3[i], f.dsub(f.dadd(f.dsub(cc[i], rr[i]), rb[i]), e2r[i]))).collect();
    let ud = uprime(&c, rr, &bxs, &cy, len);
    let rud = c.mul(rr, &ud, len);
    for i in 3..len {
        let h = f.dsub(f.dadd(rr[i], u[i]), rud[i]);
        println!("{i} {} {}", h.a, h.b);
    }
    eprintln!("done {:.1}s", t0.elapsed().as_secs_f64());
}
