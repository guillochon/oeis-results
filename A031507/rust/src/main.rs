//! A031507 / A031508: an unconditional upper bound on the rank of y^2 = x^3 + v for every v in a
//! range, by descent via 3-isogeny.
//!
//! For E_v: y^2 = x^3 + v, the 3-isogeny φ with kernel {O, (0, ±√v)} and its dual have Selmer
//! groups inside the "minus parts" of K*/K*^3 for K = Q(√v) and K = Q(√−3v) (Kummer theory for
//! Z/3 twisted by a quadratic character). Classes are unramified outside S = {p | 6v}. The minus
//! part of K(S,3) has F_3-dimension at most
//!     s(K) = h3(K) + u(K) + #{p ∈ S : p splits in K},
//! where h3 is the 3-rank of the class group and u(K) = 1 if K is real or K = Q(√−3), else 0.
//! When K = Q × Q (v a square), the group is Q(S,3), of dimension #S. Since
//!     3^r = |Im α| |Im α̂| / (|E(Q)[3]| · |E'(Q)[φ̂]/φ(E(Q)[3])|),
//! we get  rank ≤ B(v) := s(Q(√v)) + s(Q(√−3v)) − δ,  with δ = 1 iff v or −3v is a square.
//! h3 comes from an unconditional table (cubic-field counts via Hasse's theorem; tools/h3chunk.gp).
//!
//! Usage:
//!   a031507 table MAX               read data/h3/h3_*.txt, write data/h3pos.bin, data/h3neg.bin
//!   a031507 filter K N [threads]    for 1 <= k <= K and v = +k, -k: histogram of B(v); write the
//!                                   sixth-power-free k with B >= N to data/cand_{plus,minus}_N_K.txt

use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

fn data_dir() -> String {
    format!("{}/../data", env!("CARGO_MANIFEST_DIR"))
}

/// The unconditional h3 table: h3 of the fundamental discriminant ±d, indexed by d = |D|.
struct H3 {
    pos: Vec<u8>,
    neg: Vec<u8>,
}

impl H3 {
    fn load() -> Self {
        let dir = data_dir();
        let pos = fs::read(format!("{dir}/h3pos.bin")).expect("run `a031507 table` first");
        let neg = fs::read(format!("{dir}/h3neg.bin")).expect("run `a031507 table` first");
        H3 { pos, neg }
    }
    fn get(&self, d: i64) -> u32 {
        let a = d.unsigned_abs() as usize;
        let t = if d > 0 { &self.pos } else { &self.neg };
        assert!(a < t.len(), "discriminant {d} outside the h3 table");
        t[a] as u32
    }
    fn max(&self) -> u64 {
        (self.pos.len().min(self.neg.len()) - 1) as u64
    }
}

fn cmd_table(max: usize) {
    let t0 = Instant::now();
    let dir = data_dir();
    let mut pos = vec![0u8; max + 1];
    let mut neg = vec![0u8; max + 1];
    let mut entries = 0u64;
    let mut files: Vec<_> = fs::read_dir(format!("{dir}/h3"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().map_or(false, |x| x == "txt"))
        .collect();
    files.sort();
    // coverage check: the .done files must tile [1, max]
    let mut covered: Vec<(u64, u64)> = Vec::new();
    for f in &files {
        let done = fs::read_to_string(format!("{}.done", f.display())).expect("chunk without .done file");
        let get = |key: &str| -> u64 {
            let s = done.split_whitespace().find(|w| w.starts_with(key)).unwrap();
            s[key.len()..].parse().unwrap()
        };
        assert!(done.contains("bad=0"), "chunk {} has a non-power-of-3 count", f.display());
        covered.push((get("lo="), get("hi=")));
        for line in BufReader::new(fs::File::open(f).unwrap()).lines() {
            let line = line.unwrap();
            if line.is_empty() {
                continue;
            }
            let mut it = line.split(' ');
            let d: i64 = it.next().unwrap().parse().unwrap();
            let h: u8 = it.next().unwrap().parse().unwrap();
            let a = d.unsigned_abs() as usize;
            if a <= max {
                if d > 0 { pos[a] = h } else { neg[a] = h }
            }
            entries += 1;
        }
    }
    covered.sort();
    let mut next = 1u64;
    for &(lo, hi) in &covered {
        assert!(lo == next, "gap in the h3 table before {lo}");
        next = hi + 1;
    }
    assert!(next > max as u64, "h3 table covers only up to {}", next - 1);
    fs::write(format!("{dir}/h3pos.bin"), &pos).unwrap();
    fs::write(format!("{dir}/h3neg.bin"), &neg).unwrap();
    let hist = |v: &[u8]| {
        let mut h = [0u64; 8];
        for &x in v {
            h[x as usize] += 1;
        }
        h
    };
    println!(
        "table to {max}: {entries} entries; h3 histogram (D>0) {:?}, (D<0) {:?}; {:.1}s",
        hist(&pos),
        hist(&neg),
        t0.elapsed().as_secs_f64()
    );
}

/// Smallest prime factor sieve.
fn spf_sieve(n: usize) -> Vec<u32> {
    let mut spf = vec![0u32; n + 1];
    for i in 2..=n {
        if spf[i] == 0 {
            let mut j = i;
            while j <= n {
                if spf[j] == 0 {
                    spf[j] = i as u32;
                }
                j += i;
            }
        }
    }
    spf
}

/// Jacobi symbol (a/n), n odd positive.
fn jacobi(mut a: u64, mut n: u64) -> i32 {
    a %= n;
    let mut t = 1;
    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            let r = n % 8;
            if r == 3 || r == 5 {
                t = -t;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            t = -t;
        }
        a %= n;
    }
    if n == 1 { t } else { 0 }
}

/// Kronecker symbol (D/p) for a prime p.
fn kron(d: i64, p: u64) -> i32 {
    if p == 2 {
        return match d.rem_euclid(8) {
            1 | 7 => 1,
            3 | 5 => -1,
            _ => 0,
        };
    }
    jacobi(d.rem_euclid(p as i64) as u64, p)
}

/// Fundamental discriminant of Q(√(sign·m)) for squarefree m ≥ 1 (1 if sign·m = 1).
fn fund_disc(sign: i64, m: u64) -> i64 {
    let d = sign * m as i64;
    if d == 1 {
        1
    } else if d.rem_euclid(4) == 1 {
        d
    } else {
        4 * d
    }
}

/// s(K) for K = Q(√D), with S = the primes in `s_primes`.
fn s_field(d: i64, s_primes: &[u64], h3: &H3) -> u32 {
    if d == 1 {
        return s_primes.len() as u32;
    }
    let u = (d > 0 || d == -3) as u32;
    let split = s_primes.iter().filter(|&&p| kron(d, p) == 1).count() as u32;
    h3.get(d) + u + split
}

/// Rank bound B(v) for v = sign·k, or None if k is not sixth-power free.
fn bound(sign: i64, k: u64, spf: &[u32], h3: &H3, s_primes: &mut Vec<u64>) -> Option<u32> {
    parts(sign, k, spf, h3, s_primes).map(|(a, b, d)| a + b - d)
}

/// (s(Q(√−3v)), s(Q(√v)), δ): the ambient bounds for Sel^φ(E) and Sel^φ̂(E'), and δ.
fn parts(sign: i64, k: u64, spf: &[u32], h3: &H3, s_primes: &mut Vec<u64>) -> Option<(u32, u32, u32)> {
    s_primes.clear();
    s_primes.push(2);
    s_primes.push(3);
    let mut m = 1u64; // squarefree kernel of k (primes with odd exponent)
    let mut x = k;
    while x > 1 {
        let p = spf[x as usize] as u64;
        let mut e = 0;
        while x % p == 0 {
            x /= p;
            e += 1;
        }
        if e >= 6 {
            return None;
        }
        if e % 2 == 1 {
            m *= p;
        }
        if p > 3 {
            s_primes.push(p);
        }
    }
    let d1 = fund_disc(sign, m);
    // squarefree part of 3m, with sign flipped (field Q(√−3v))
    let m3 = if m % 3 == 0 { m / 3 } else { 3 * m };
    let d2 = fund_disc(-sign, m3);
    let delta = (d1 == 1 || d2 == 1) as u32;
    Some((s_field(d2, s_primes, h3), s_field(d1, s_primes, h3), delta))
}

fn cmd_filter(kmax: u64, n: u32, threads: usize) {
    let t0 = Instant::now();
    let h3 = H3::load();
    assert!(12 * kmax <= h3.max(), "h3 table too small: need 12K = {}", 12 * kmax);
    let spf = spf_sieve(kmax as usize);
    eprintln!("setup {:.1}s", t0.elapsed().as_secs_f64());
    let next = AtomicU64::new(1);
    const CHUNK: u64 = 1 << 16;
    let results: Vec<([[u64; 32]; 2], [Vec<(u64, u32, u32, u32)>; 2])> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads)
            .map(|_| {
                sc.spawn(|| {
                    let mut hist = [[0u64; 32]; 2];
                    let mut cand: [Vec<(u64, u32, u32, u32)>; 2] = [Vec::new(), Vec::new()];
                    let mut sp = Vec::with_capacity(16);
                    loop {
                        let lo = next.fetch_add(CHUNK, Ordering::Relaxed);
                        if lo > kmax {
                            break;
                        }
                        let hi = (lo + CHUNK - 1).min(kmax);
                        for k in lo..=hi {
                            for (si, sign) in [(0usize, 1i64), (1, -1)] {
                                if let Some((sa, sb, d)) = parts(sign, k, &spf, &h3, &mut sp) {
                                    let b = sa + sb - d;
                                    hist[si][b as usize] += 1;
                                    if b >= n {
                                        cand[si].push((k, sa, sb, d));
                                    }
                                }
                            }
                        }
                    }
                    (hist, cand)
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut hist = [[0u64; 32]; 2];
    let mut cand: [Vec<(u64, u32, u32, u32)>; 2] = [Vec::new(), Vec::new()];
    for (h, c) in results {
        for s in 0..2 {
            for b in 0..32 {
                hist[s][b] += h[s][b];
            }
            cand[s].extend(c[s].iter());
        }
    }
    let dir = data_dir();
    for (s, name) in [(0, "plus"), (1, "minus")] {
        cand[s].sort();
        let path = format!("{dir}/cand_{name}_{n}_{kmax}.txt");
        let mut f = fs::File::create(&path).unwrap();
        // columns: k, s_a = s(Q(√−3v)) (bounds Sel^φ), s_b = s(Q(√v)) (bounds Sel^φ̂), δ
        for (k, sa, sb, d) in &cand[s] {
            writeln!(f, "{k} {sa} {sb} {d}").unwrap();
        }
        let top = hist[s].iter().rposition(|&x| x > 0).unwrap_or(0);
        println!(
            "{name}: k <= {kmax}, B histogram {:?}; {} with B >= {n} -> {path}",
            &hist[s][..=top],
            cand[s].len()
        );
    }
    println!("time {:.1}s", t0.elapsed().as_secs_f64());
}

/// Write "k B(+k) B(-k)" for 1 <= k <= K (B = -1 if k is not sixth-power free), for cross-checks.
fn cmd_dump(kmax: u64) {
    let h3 = H3::load();
    let spf = spf_sieve(kmax as usize);
    let mut sp = Vec::new();
    let mut out = String::new();
    for k in 1..=kmax {
        let b = |s: i64, sp: &mut Vec<u64>| bound(s, k, &spf, &h3, sp).map_or(-1, |x| x as i64);
        let (bp, bm) = (b(1, &mut sp), b(-1, &mut sp));
        out.push_str(&format!("{k} {bp} {bm}\n"));
    }
    let path = format!("{}/dump_{kmax}.txt", data_dir());
    fs::write(&path, out).unwrap();
    println!("wrote {path}");
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
    let num = |i: usize, d: u64| args.get(i).map_or(d, |s| s.parse().expect("number"));
    match args.get(1).map(String::as_str) {
        Some("table") => cmd_table(num(2, 683_000_000) as usize),
        Some("dump") => cmd_dump(num(2, 100000)),
        Some("filter") => cmd_filter(num(2, 0), num(3, 7) as u32, num(4, 10) as usize),
        _ => eprintln!("usage: a031507 table MAX | filter K N [threads]"),
    }
}
