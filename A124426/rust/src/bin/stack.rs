//! Stacking matryoshkas with distinct sizes: the slot recurrence of stack_fast.py on a dense array.
//!
//! State (s, f, e, t) after k dolls: s free heads at tower tops, f towers ending in an empty open bottom,
//! e outer closed dolls (each holds one chain), t outer tops (each holds one chain of tops).
//! Reachable states satisfy e + f <= k, e + t <= k, s <= e + t.
//!
//! One doll maps the state vector A to D(A) + T(B(A)), where D, B, T place a closed doll, a bottom and a
//! top. Each is written as a gather and applied in place, sweeping the coordinate it lowers in descending
//! order, so only two arrays are needed.
//!
//!   stack N --log2                 log2 a(n) for n <= N (f64, rescaled each step)
//!   stack N --primes FIRST COUNT   a(n) mod each of the primes below 2^31 with index FIRST..FIRST+COUNT
//!                                  (prints "n p residue"; reassemble with the Chinese remainder theorem)

use matryoshka::{disable_power_throttling, primes_below, threads};
use std::time::Instant;

struct Layout {
    n: usize,
    base: Vec<usize>, // (e, f, t) -> offset of s = 0; usize::MAX if invalid
    slab: Vec<usize>, // start of slab e; slab[n + 1] = total length
}

impl Layout {
    fn new(n: usize) -> Self {
        let m = n + 1;
        let mut base = vec![usize::MAX; m * m * m];
        let mut slab = vec![0; m + 1];
        let mut off = 0;
        for e in 0..=n {
            slab[e] = off;
            for f in 0..=n - e {
                for t in 0..=n - e {
                    base[(e * m + f) * m + t] = off;
                    off += (e + t).min(n) + 1;
                }
            }
        }
        slab[m] = off;
        Layout { n, base, slab }
    }
    #[inline]
    fn smax(&self, e: usize, t: usize) -> usize {
        (e + t).min(self.n)
    }
    #[inline]
    fn at(&self, e: usize, f: usize, t: usize) -> usize {
        let m = self.n + 1;
        self.base[(e * m + f) * m + t]
    }
}

trait Val: Copy + Send + Sync + Default {
    fn add(self, o: Self, p: u64) -> Self;
    fn mul(self, k: usize, p: u64) -> Self;
}

impl Val for u32 {
    #[inline]
    fn add(self, o: Self, p: u64) -> Self {
        let s = self as u64 + o as u64;
        (if s >= p { s - p } else { s }) as u32
    }
    #[inline]
    fn mul(self, k: usize, p: u64) -> Self {
        ((self as u64 * k as u64) % p) as u32
    }
}

impl Val for f64 {
    #[inline]
    fn add(self, o: Self, _p: u64) -> Self {
        self + o
    }
    #[inline]
    fn mul(self, k: usize, _p: u64) -> Self {
        self * k as f64
    }
}

/// Split `arr` into the slabs e = 0..=k (each a mutable slice), using the layout's slab starts.
fn slabs_mut<'a, T>(arr: &'a mut [T], lay: &Layout, k: usize) -> Vec<&'a mut [T]> {
    let mut out = Vec::with_capacity(k + 1);
    let mut rest = &mut arr[..lay.slab[k + 1]];
    for e in 0..=k {
        let len = lay.slab[e + 1] - lay.slab[e];
        let (head, tail) = rest.split_at_mut(len);
        out.push(head);
        rest = tail;
    }
    out
}

/// Run `work(e, slab)` for every slab e = 0..=k, spread over the thread pool.
fn par_slabs<T: Send, F: Fn(usize, &mut [T]) + Sync>(slabs: Vec<&mut [T]>, nthreads: usize, work: F) {
    let mut buckets: Vec<Vec<(usize, &mut [T])>> = (0..nthreads).map(|_| Vec::new()).collect();
    // largest slabs first, round robin, for balance
    let mut items: Vec<(usize, &mut [T])> = slabs.into_iter().enumerate().collect();
    items.sort_by_key(|(_, s)| std::cmp::Reverse(s.len()));
    for (i, it) in items.into_iter().enumerate() {
        buckets[i % nthreads].push(it);
    }
    std::thread::scope(|sc| {
        for b in buckets {
            let w = &work;
            sc.spawn(move || {
                for (e, s) in b {
                    w(e, s);
                }
            });
        }
    });
}

fn step<T: Val>(lay: &Layout, k: usize, a: &mut [T], tmp: &mut [T], p: u64, nthreads: usize) {
    let z = T::default();
    // 1. tmp = a on the region of k dolls
    {
        let src: &[T] = a;
        par_slabs(slabs_mut(tmp, lay, k), nthreads, |e, sl| {
            let s0 = lay.slab[e];
            for f in 0..=k - e {
                let lo = lay.at(e, f, 0);
                let hi = lay.at(e, f, k - e) + lay.smax(e, k - e).min(k) + 1;
                sl[lo - s0..hi - s0].copy_from_slice(&src[lo..hi]);
            }
        });
    }
    // 2. B in place on tmp (f descending): new = (f+e)*old[f] + old[s,f-1] + (s+1)*old[s+1,f-1]
    par_slabs(slabs_mut(tmp, lay, k), nthreads, |e, sl| {
        let s0 = lay.slab[e];
        for f in (0..=k - e).rev() {
            for t in 0..=k - e {
                let sm = lay.smax(e, t).min(k);
                let cur = lay.at(e, f, t) - s0;
                let prv = if f > 0 { lay.at(e, f - 1, t) - s0 } else { usize::MAX };
                for s in 0..=sm {
                    let mut v = sl[cur + s].mul(f + e, p);
                    if f > 0 {
                        v = v.add(sl[prv + s], p);
                        if s < sm {
                            v = v.add(sl[prv + s + 1].mul(s + 1, p), p);
                        }
                    }
                    sl[cur + s] = v;
                }
            }
        }
    });
    // 3. T in place on tmp (t descending): new = t*old[t] + old[s-1,t-1] + s*old[s,t-1]
    par_slabs(slabs_mut(tmp, lay, k), nthreads, |e, sl| {
        let s0 = lay.slab[e];
        for f in 0..=k - e {
            for t in (0..=k - e).rev() {
                let sm = lay.smax(e, t).min(k);
                let cur = lay.at(e, f, t) - s0;
                if t == 0 {
                    for s in 0..=sm {
                        sl[cur + s] = z;
                    }
                    continue;
                }
                let prv = lay.at(e, f, t - 1) - s0;
                let smp = lay.smax(e, t - 1).min(k);
                for s in (0..=sm).rev() {
                    let mut v = sl[cur + s].mul(t, p);
                    if s >= 1 && s - 1 <= smp {
                        v = v.add(sl[prv + s - 1], p);
                    }
                    if s <= smp {
                        v = v.add(sl[prv + s].mul(s, p), p);
                    }
                    sl[cur + s] = v;
                }
            }
        }
    });
    // 4. D in place on a (e descending; slab e reads slab e-1):
    //    new = e*old[e] + old[s-1,e-1] + s*old[s,e-1] + (f+1)*old[s-1,f+1,e-1]
    for e in (0..=k).rev() {
        let (lower, upper) = a.split_at_mut(lay.slab[e]);
        let slab = &mut upper[..lay.slab[e + 1] - lay.slab[e]];
        let s0 = lay.slab[e];
        let lower: &[T] = lower;
        // split the slab by f
        let mut blocks: Vec<(usize, &mut [T])> = Vec::new();
        let mut rest = slab;
        let mut pos = s0;
        for f in 0..=lay.n - e {
            let end = if f < lay.n - e { lay.at(e, f + 1, 0) } else { lay.slab[e + 1] };
            let (h, tl) = rest.split_at_mut(end - pos);
            if f <= k - e {
                blocks.push((f, h));
            }
            rest = tl;
            pos = end;
        }
        let nt = nthreads.min(blocks.len()).max(1);
        let mut buckets: Vec<Vec<(usize, &mut [T])>> = (0..nt).map(|_| Vec::new()).collect();
        for (i, b) in blocks.into_iter().enumerate() {
            buckets[i % nt].push(b);
        }
        std::thread::scope(|sc| {
            for bk in buckets {
                sc.spawn(move || {
                    for (f, blk) in bk {
                        let b0 = lay.at(e, f, 0);
                        for t in 0..=k - e {
                            let sm = lay.smax(e, t).min(k);
                            let cur = lay.at(e, f, t) - b0;
                            for s in 0..=sm {
                                let mut v = blk[cur + s].mul(e, p);
                                if e >= 1 {
                                    let smp = lay.smax(e - 1, t);
                                    let q = lay.at(e - 1, f, t);
                                    if s >= 1 && s - 1 <= smp {
                                        v = v.add(lower[q + s - 1], p);
                                    }
                                    if s <= smp {
                                        v = v.add(lower[q + s].mul(s, p), p);
                                    }
                                    if s >= 1 && s - 1 <= smp && f + 1 <= lay.n - (e - 1) {
                                        let q2 = lay.at(e - 1, f + 1, t);
                                        v = v.add(lower[q2 + s - 1].mul(f + 1, p), p);
                                    }
                                }
                                blk[cur + s] = v;
                            }
                        }
                    }
                });
            }
        });
        let _ = s0;
    }
    // 5. a += tmp on the region
    {
        let src: &[T] = tmp;
        par_slabs(slabs_mut(a, lay, k), nthreads, |e, sl| {
            let s0 = lay.slab[e];
            for f in 0..=k - e {
                let lo = lay.at(e, f, 0);
                let hi = lay.at(e, f, k - e) + lay.smax(e, k - e).min(k) + 1;
                for i in lo..hi {
                    sl[i - s0] = sl[i - s0].add(src[i], p);
                }
            }
        });
    }
}

fn total<T: Val>(lay: &Layout, k: usize, a: &[T], p: u64) -> T {
    let mut acc = T::default();
    for e in 0..=k {
        for f in 0..=k - e {
            let lo = lay.at(e, f, 0);
            let hi = lay.at(e, f, k - e) + lay.smax(e, k - e).min(k) + 1;
            for x in &a[lo..hi] {
                acc = acc.add(*x, p);
            }
        }
    }
    acc
}

fn main() {
    disable_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let n: usize = args[0].parse().expect("N");
    let lay = Layout::new(n);
    let len = lay.slab[n + 1];
    let nthreads = threads();
    eprintln!("N = {n}: {len} states per array ({:.2} GB per array of u32)", len as f64 * 4.0 / 1e9);
    if args.get(1).map(String::as_str) == Some("--log2") {
        let mut a = vec![0f64; len];
        let mut tmp = vec![0f64; len];
        a[lay.at(0, 0, 0)] = 1.0;
        let mut logscale = 0f64;
        for k in 1..=n {
            step(&lay, k, &mut a, &mut tmp, 0, nthreads);
            let tot = total(&lay, k, &a, 0);
            println!("{k} {:.6}", logscale + tot.log2());
            // rescale so values stay in range
            let lg = tot.log2().floor();
            let sc = (2f64).powf(-lg);
            for x in a.iter_mut() {
                *x *= sc;
            }
            logscale += lg;
        }
        return;
    }
    assert_eq!(args.get(1).map(String::as_str), Some("--primes"));
    let first: usize = args[2].parse().unwrap();
    let count: usize = args[3].parse().unwrap();
    let primes = primes_below(1u64 << 31, first + count);
    let mut a = vec![0u32; len];
    let mut tmp = vec![0u32; len];
    for &p in &primes[first..] {
        a.iter_mut().for_each(|x| *x = 0);
        a[lay.at(0, 0, 0)] = 1;
        let t0 = Instant::now();
        for k in 1..=n {
            let ts = Instant::now();
            step(&lay, k, &mut a, &mut tmp, p, nthreads);
            println!("{k} {p} {}", total(&lay, k, &a, p));
            eprintln!("step {k} p {p} {:.3}s", ts.elapsed().as_secs_f64());
        }
        eprintln!("prime {p} done in {:.1}s", t0.elapsed().as_secs_f64());
    }
}
