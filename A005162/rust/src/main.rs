//! A005162 (DASASMs) modulo 31-bit primes: the column-sweep six-vertex transfer matrix of `dasasm.py`,
//! in place on one u32 array, multithreaded.
//!
//! Layout: bit k (0-based) = running sum of stored open row k+1; while a column is processed, bit R (R = number of
//! stored rows) is the vertical running sum v. For even n the middle row n/2 opens at column n/2 and closes at
//! column n/2+1 without being touched in between, so it is not stored (its two values give a factor 2).
//! Peak memory: 2^ceil(n/2) words.
//!
//! Usage: dasasm <n> <k0> <k1> [--threads T]   residues of a(n) modulo primes number k0..k1-1 (largest first
//! below 2^31). Prints "n p residue seconds" per prime.

use std::time::Instant;

// Storage word. The u16 build (`--features w16`) halves memory and uses primes below 2^16; arithmetic is in u32.
#[cfg(not(feature = "w16"))]
type W = u32;
#[cfg(not(feature = "w16"))]
const PBITS: u32 = 31;
#[cfg(feature = "w16")]
type W = u16;
#[cfg(feature = "w16")]
const PBITS: u32 = 16;

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut b: u64, mut e: u64| {
        let mut r = 1u64;
        b %= n;
        while e > 0 {
            if e & 1 == 1 {
                r = mulmod(r, b);
            }
            b = mulmod(b, b);
            e >>= 1;
        }
        r
    };
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if a % n == 0 {
            continue;
        }
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn primes(k: usize) -> Vec<u32> {
    let mut out = Vec::new();
    let mut q: u64 = (1 << PBITS) - 1;
    while out.len() < k {
        if is_prime(q) {
            out.push(q as u32);
        }
        q -= 2;
    }
    out
}

#[derive(Clone, Copy)]
struct Ptr(*mut W);
unsafe impl Send for Ptr {}
unsafe impl Sync for Ptr {}
impl Ptr {
    #[inline(always)]
    fn at(self, i: usize) -> *mut W {
        unsafe { self.0.add(i) }
    }
}

/// Run f(i) for i in 0..len over `threads` contiguous ranges.
fn par_for<F: Fn(usize) + Sync>(len: usize, threads: usize, f: F) {
    if len < (1 << 16) || threads == 1 {
        (0..len).for_each(&f);
        return;
    }
    let chunk = len.div_ceil(threads);
    std::thread::scope(|s| {
        for t in 0..threads {
            let (lo, hi) = (t * chunk, ((t + 1) * chunk).min(len));
            let f = &f;
            s.spawn(move || (lo..hi).for_each(f));
        }
    });
}

#[inline(always)]
fn add(a: W, b: W, p: u32) -> W {
    let s = a as u32 + b as u32; // a, b < 2^31 so no overflow
    (if s >= p { s - p } else { s }) as W
}

fn count_mod(n: usize, p: u32, a: &mut [W], threads: usize) -> u32 {
    let skip = if n % 2 == 0 { n / 2 } else { 0 }; // unstored middle row (even n)
    let mut factor: u64 = 1;
    let mut r = 0usize; // stored open rows
    a[0] = 1;
    for c in 1..=n {
        let h = c.min(n + 1 - c);
        let half = 1usize << r;
        // v = 0: clear the upper half
        let ptr = Ptr(a.as_mut_ptr());
        par_for(half, threads, |i| unsafe { *ptr.at(half + i) = 0 });
        // interior rows 1..h-1 are stored bits 0..h-2 (rows below the middle are always stored)
        for row in 0..h - 1 {
            let rb = 1usize << row;
            let low = rb - 1;
            par_for(half / 2, threads, |q| unsafe {
                let i = ((q & !low) << 1) | (q & low); // insert a 0 at bit `row`
                let j = half + (i | rb);
                let s = add(*ptr.at(i), *ptr.at(j), p);
                *ptr.at(i) = s;
                *ptr.at(j) = s;
            });
        }
        if c < n + 1 - c {
            // diagonal: drop v, open row c (free bit)
            if c == skip {
                par_for(half, threads, |i| unsafe { *ptr.at(i) = add(*ptr.at(i), *ptr.at(half + i), p) });
                factor = factor * 2 % p as u64;
            } else {
                par_for(half, threads, |i| unsafe {
                    let s = add(*ptr.at(i), *ptr.at(half + i), p);
                    *ptr.at(i) = s;
                    *ptr.at(half + i) = s;
                });
                r += 1;
            }
        } else if c > n + 1 - c && n + 1 - c != skip {
            // antidiagonal: close the top stored row (bit r-1) and drop v
            let q4 = half / 2;
            par_for(q4, threads, |i| unsafe {
                let s = (*ptr.at(i) as u64
                    + *ptr.at(i + q4) as u64
                    + *ptr.at(half + i) as u64
                    + *ptr.at(half + i + q4) as u64)
                    % p as u64;
                *ptr.at(i) = s as W;
            });
            r -= 1;
        } else {
            // centre (odd n) or closing the unstored middle row: drop v
            par_for(half, threads, |i| unsafe { *ptr.at(i) = add(*ptr.at(i), *ptr.at(half + i), p) });
        }
    }
    assert_eq!(r, 0);
    (a[0] as u64 * factor % p as u64) as u32
}

#[cfg(windows)]
fn disable_power_throttling() {
    // Without this, background jobs on this PC are scheduled on the E-cores only (EcoQoS).
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
    let ok = unsafe { SetProcessInformation(GetCurrentProcess(), 4, &st as *const _ as *const _, 12) };
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
    let threads: usize = opt("--threads").map(|s| s.parse().unwrap()).unwrap_or(10);
    let n: usize = args[0].parse().expect("usage: dasasm <n> <k0> <k1> [--threads T]");
    let k0: usize = args[1].parse().unwrap();
    let k1: usize = args[2].parse().unwrap();
    let ps = primes(k1);
    let peak = 1usize << n.div_ceil(2);
    let mut a: Vec<W> = vec![0; peak.max(2)];
    for &p in &ps[k0..k1] {
        let t0 = Instant::now();
        a.iter_mut().for_each(|x| *x = 0);
        let res = count_mod(n, p, &mut a, threads);
        println!("{} {} {} {:.1}", n, p, res, t0.elapsed().as_secs_f64());
    }
}
