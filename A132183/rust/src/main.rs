//! OEIS A132183: regular Boolean functions of n variables = order ideals of the binary
//! majorization lattice M(n) (subsets of [n]; A <= B iff B is reached from A by adding
//! elements and moving elements up).
//!
//! Split n = m + k. A subset is (A, P) with A in M(m) and a pattern P in M(k) on the top k
//! elements. An ideal of M(n) is a family (I_P) of ideals of M(m), one per pattern, with
//!   I_P >= I_P'          when P' is P with one element moved up,
//!   I_P >= sigma(I_P')   when P' is P plus the lowest pattern element,
//! where sigma(I) = {A : A minus its largest element is in I}.
//! With L = J(M(m)), sums over one I_P are zeta transforms on the distributive lattice L,
//! done in |Hasse edges of L| additions each.
//!
//!   regular-ideals k2 m [threads]         a(m+2)   (chain of 4 patterns; vectors over L)
//!   regular-ideals k4 m [threads] [prime] a(m+4)   (16 patterns; tables on comparable pairs)
//!
//! Both modes also print the number of self-dual regular functions (A109456) as a check.
//! k4 prints a(m+4) mod one 32-bit prime per run (prime index 0..3); combine with crt.py.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

type Bs = [u64; 4];
const PRIMES: [u64; 4] = [4294967291, 4294967279, 4294967231, 4294967197];

#[inline]
fn bit(b: &Bs, i: usize) -> bool {
    b[i >> 6] >> (i & 63) & 1 == 1
}
#[inline]
fn set(b: &mut Bs, i: usize) {
    b[i >> 6] |= 1u64 << (i & 63);
}
#[inline]
fn clr(b: &mut Bs, i: usize) {
    b[i >> 6] &= !(1u64 << (i & 63));
}
#[inline]
fn subset(a: &Bs, b: &Bs) -> bool {
    (a[0] & !b[0]) | (a[1] & !b[1]) | (a[2] & !b[2]) | (a[3] & !b[3]) == 0
}
#[inline]
fn meets(a: &Bs, b: &Bs) -> bool {
    (a[0] & b[0]) | (a[1] & b[1]) | (a[2] & b[2]) | (a[3] & b[3]) != 0
}

fn leq(a: u32, b: u32) -> bool {
    // dominance: i-th largest element of a <= i-th largest of b, and |a| <= |b|
    let (mut x, mut y) = (a, b);
    while x != 0 {
        if y == 0 {
            return false;
        }
        let (hx, hy) = (31 - x.leading_zeros(), 31 - y.leading_zeros());
        if hx > hy {
            return false;
        }
        x &= !(1 << hx);
        y &= !(1 << hy);
    }
    true
}

struct Lat {
    n: usize,          // 2^m elements of M(m), in a linear extension
    ideals: Vec<Bs>,   // sorted
    off: Vec<usize>,   // Hasse edges of L grouped by removed element p: off[p]..off[p+1]
    hi: Vec<u32>,
    lo: Vec<u32>,
    sig: Vec<u32>,
    dual: Vec<u32>,
}

impl Lat {
    fn new(m: usize) -> Lat {
        let n = 1usize << m;
        let rank = |a: u32| (0..m).filter(|i| a >> i & 1 == 1).map(|i| i + 1).sum::<usize>();
        let mut elems: Vec<u32> = (0..n as u32).collect();
        elems.sort_by_key(|&a| (rank(a), a));
        let mut pos = vec![0usize; n];
        for (i, &a) in elems.iter().enumerate() {
            pos[a as usize] = i;
        }
        let mut below = vec![[0u64; 4]; n];
        let mut above = vec![[0u64; 4]; n];
        for i in 0..n {
            for j in 0..i {
                if leq(elems[j], elems[i]) {
                    set(&mut below[i], j);
                    set(&mut above[j], i);
                }
            }
        }
        // ideals by DFS over the linear extension
        let mut ideals: Vec<Bs> = Vec::new();
        fn rec(i: usize, n: usize, cur: &mut Bs, below: &[Bs], out: &mut Vec<Bs>) {
            if i == n {
                out.push(*cur);
                return;
            }
            rec(i + 1, n, cur, below, out);
            if subset(&below[i], cur) {
                set(cur, i);
                rec(i + 1, n, cur, below, out);
                clr(cur, i);
            }
        }
        let mut cur = [0u64; 4];
        rec(0, n, &mut cur, &below, &mut ideals);
        ideals.sort_unstable();
        let nn = ideals.len();
        let find = |b: &Bs| ideals.binary_search(b).expect("ideal not found") as u32;
        // Hasse edges
        let mut cnt = vec![0usize; n + 1];
        for id in &ideals {
            for p in 0..n {
                if bit(id, p) && !meets(&above[p], id) {
                    cnt[p + 1] += 1;
                }
            }
        }
        let mut off = vec![0usize; n + 1];
        for p in 0..n {
            off[p + 1] = off[p] + cnt[p + 1];
        }
        let ne = off[n];
        let mut hi = vec![0u32; ne];
        let mut lo = vec![0u32; ne];
        let mut fill = off.clone();
        for (t, id) in ideals.iter().enumerate() {
            for p in 0..n {
                if bit(id, p) && !meets(&above[p], id) {
                    let mut x = *id;
                    clr(&mut x, p);
                    hi[fill[p]] = t as u32;
                    lo[fill[p]] = find(&x);
                    fill[p] += 1;
                }
            }
        }
        // sigma and dual
        let dm: Vec<usize> = elems
            .iter()
            .map(|&a| pos[if a == 0 { 0 } else { (a & !(1 << (31 - a.leading_zeros()))) as usize }])
            .collect();
        let comp: Vec<usize> = elems.iter().map(|&a| pos[(n as u32 - 1 - a) as usize]).collect();
        let mut sig = vec![0u32; nn];
        let mut dual = vec![0u32; nn];
        for (t, id) in ideals.iter().enumerate() {
            let (mut s, mut d) = ([0u64; 4], [0u64; 4]);
            for i in 0..n {
                if bit(id, dm[i]) {
                    set(&mut s, i);
                }
                if !bit(id, i) {
                    set(&mut d, comp[i]);
                }
            }
            sig[t] = find(&s);
            dual[t] = find(&d);
        }
        Lat { n, ideals, off, hi, lo, sig, dual }
    }

    /// g(K) = sum over I >= K of f(I)
    fn up<T: Copy + std::ops::AddAssign>(&self, f: &mut [T]) {
        for p in (0..self.n).rev() {
            for e in self.off[p]..self.off[p + 1] {
                let x = f[self.hi[e] as usize];
                f[self.lo[e] as usize] += x;
            }
        }
    }
    /// g(J) = sum over I <= J of f(I)
    fn down<T: Copy + std::ops::AddAssign>(&self, f: &mut [T]) {
        for p in 0..self.n {
            for e in self.off[p]..self.off[p + 1] {
                let x = f[self.lo[e] as usize];
                f[self.hi[e] as usize] += x;
            }
        }
    }
}

fn k2(l: &Lat, m: usize) {
    let nn = l.ideals.len();
    let mut c = vec![1u64; nn];
    l.up(&mut c);
    let u: Vec<u64> = (0..nn).map(|j| c[l.sig[j] as usize]).collect();
    drop(c);
    let mut w = vec![0u64; nn];
    for j in 0..nn {
        w[l.sig[j] as usize] += 1;
    }
    l.down(&mut w);
    for j in 0..nn {
        assert_eq!(w[j], u[l.dual[j] as usize]);
    }
    l.down(&mut w);
    let mut total = 0u128;
    let mut sd = 0u128;
    for j in 0..nn {
        total += u[j] as u128 * w[j] as u128;
        if subset(&l.ideals[l.dual[j] as usize], &l.ideals[j]) {
            sd += u[j] as u128;
        }
    }
    println!("k2 m={} a({}) = {}", m, m + 2, total);
    println!("k2 m={} selfdual({}) = {}", m, m + 2, sd);
}

/// Row d of the tables: P(d, e) and R(d, g) = sum over e in [g, d] of P(d, e).
struct Work {
    a: Vec<u128>,
    b: Vec<u128>,
}

fn row(l: &Lat, v: &[u128], d: usize, w: &mut Work) -> u128 {
    let nn = l.ideals.len();
    let idd = &l.ideals[d];
    let a = &mut w.a;
    let b = &mut w.b;
    // n_d(I_2) = #{I_12 >= d : sigma(I_12) <= I_2}
    a.iter_mut().for_each(|x| *x = 0);
    for i in 0..nn {
        if subset(idd, &l.ideals[i]) {
            a[l.sig[i] as usize] += 1;
        }
    }
    l.down(a);
    // T2(I_3, d) = sum over I_2 >= I_3 of v(I_2) n_d(I_2), kept only for I_3 >= sigma(d)
    for i in 0..nn {
        a[i] *= v[i];
    }
    l.up(a);
    let sd_id = &l.ideals[l.sig[d] as usize];
    for i in 0..nn {
        if !subset(sd_id, &l.ideals[i]) {
            a[i] = 0;
        }
    }
    l.up(a); // T3(I_4, d)
    l.up(a); // sum over I_4 >= K
    // P(d, e) = that at K = sigma(e), for e <= d
    let dd = &l.ideals[l.dual[d] as usize]; // e* <= d  iff  d* <= e
    let mut sd = 0u128;
    for e in 0..nn {
        if subset(&l.ideals[e], idd) {
            let x = a[l.sig[e] as usize];
            b[e] = x;
            if subset(dd, &l.ideals[e]) {
                sd += x;
            }
        } else {
            b[e] = 0;
        }
    }
    l.up(b); // R(d, g)
    sd
}

fn k4(l: &Lat, m: usize, threads: usize, pi: usize) {
    let nn = l.ideals.len();
    let p = PRIMES[pi];
    let t0 = Instant::now();
    let mut c = vec![1u128; nn];
    l.up(&mut c);
    let mut v: Vec<u128> = (0..nn).map(|j| c[l.sig[j] as usize]).collect();
    l.up(&mut v);
    // column sizes: number of d >= g
    let mut col = vec![0usize; nn + 1];
    for g in 0..nn {
        col[g + 1] = col[g] + c[g] as usize;
    }
    let total_pairs = col[nn];
    eprintln!("comparable pairs: {}", total_pairs);
    let mut t = vec![0u32; total_pairs]; // t[col[g] + k] = R(d_k, g) mod p, d_k = k-th ideal >= g
    let mut cursor: Vec<usize> = col[..nn].to_vec();
    // phase 1: rows in chunks, scattered in row order
    let chunk = 64 * threads;
    let mut sd_total = 0u128;
    let mut start = 0;
    while start < nn {
        let end = (start + chunk).min(nn);
        let next = AtomicUsize::new(start);
        let mut rows: Vec<(usize, Vec<u32>, u128)> = Vec::new();
        std::thread::scope(|s| {
            let hs: Vec<_> = (0..threads)
                .map(|_| {
                    s.spawn(|| {
                        let mut w = Work { a: vec![0; nn], b: vec![0; nn] };
                        let mut out = Vec::new();
                        loop {
                            let d = next.fetch_add(1, Ordering::Relaxed);
                            if d >= end {
                                break;
                            }
                            let sd = row(l, &v, d, &mut w);
                            let idd = &l.ideals[d];
                            let vals: Vec<u32> = (0..nn)
                                .filter(|&g| subset(&l.ideals[g], idd))
                                .map(|g| (w.b[g] % p as u128) as u32)
                                .collect();
                            out.push((d, vals, sd));
                        }
                        out
                    })
                })
                .collect();
            for h in hs {
                rows.extend(h.join().unwrap());
            }
        });
        rows.sort_by_key(|r| r.0);
        for (d, vals, sd) in rows {
            sd_total += sd;
            let idd = &l.ideals[d];
            let mut k = 0;
            for g in 0..nn {
                if subset(&l.ideals[g], idd) {
                    t[cursor[g]] = vals[k];
                    cursor[g] += 1;
                    k += 1;
                }
            }
            assert_eq!(k, vals.len());
        }
        start = end;
    }
    for g in 0..nn {
        assert_eq!(cursor[g], col[g + 1]);
    }
    eprintln!("phase 1 done in {:.1} s", t0.elapsed().as_secs_f64());
    println!("k4 m={} selfdual({}) = {}", m, m + 4, sd_total);
    // phase 2: sum over d >= g of R(d, g) R(g*, d*); column d* of t lists R(x, d*), x >= d*
    let next = AtomicUsize::new(0);
    let mut acc = 0u64;
    let tref = &t;
    let colref = &col;
    let vref = &v;
    std::thread::scope(|s| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut w = Work { a: vec![0; nn], b: vec![0; nn] };
                    let mut acc = 0u64;
                    loop {
                        let d = next.fetch_add(1, Ordering::Relaxed);
                        if d >= nn {
                            break;
                        }
                        row(l, vref, d, &mut w);
                        let y = l.dual[d] as usize;
                        let idy = &l.ideals[y];
                        let mut k = colref[y];
                        for x in 0..nn {
                            if subset(idy, &l.ideals[x]) {
                                let r = (w.b[l.dual[x] as usize] % p as u128) as u64;
                                acc = (acc + r * tref[k] as u64 % p) % p;
                                k += 1;
                            }
                        }
                        assert_eq!(k, colref[y + 1]);
                    }
                    acc
                })
            })
            .collect();
        for h in hs {
            acc = (acc + h.join().unwrap()) % p;
        }
    });
    println!("k4 m={} a({}) mod {} = {}", m, m + 4, p, acc);
    eprintln!("done in {:.1} s", t0.elapsed().as_secs_f64());
}

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
    let st = ProcessPowerThrottlingState { version: 1, control_mask: 1, state_mask: 0 };
    unsafe {
        SetProcessInformation(
            GetCurrentProcess(),
            4,
            &st as *const _ as *const std::ffi::c_void,
            std::mem::size_of::<ProcessPowerThrottlingState>() as u32,
        );
    }
}
#[cfg(not(windows))]
fn disable_power_throttling() {}

fn main() {
    disable_power_throttling();
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: regular-ideals k2|k4 m [threads] [prime index 0..3]");
        std::process::exit(2);
    }
    let m: usize = args[2].parse().unwrap();
    assert!(m <= 8);
    let threads: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(6);
    let pi: usize = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(0);
    let t0 = Instant::now();
    let l = Lat::new(m);
    eprintln!(
        "m={}: {} ideals, {} Hasse edges, built in {:.1} s",
        m,
        l.ideals.len(),
        l.hi.len(),
        t0.elapsed().as_secs_f64()
    );
    match args[1].as_str() {
        "k2" => k2(&l, m),
        "k4" => k4(&l, m, threads, pi),
        _ => panic!("unknown mode"),
    }
}
