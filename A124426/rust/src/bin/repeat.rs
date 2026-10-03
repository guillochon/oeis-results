//! Repeated sizes (see repeat_fast.py for the derivation): Burnside over permutations of equal pieces.
//! Fix(g) = prod_L Lab_L, where Lab_L is the labelled slot recurrence on the L-cycles of g, with every
//! non-table placement weighted L. One DP over the parts of the composition (sizes, largest first) sums
//! over compositions, closed/open splits and cycle types (weight 1/(z_D z_B z_T)).
//! Values are kept modulo P primes near 2^61 and printed per n as "n r_1 ... r_P".
//!
//!   repeat nest|stack N

use matryoshka::{disable_power_throttling, mulmod, powmod, primes_below, threads};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

const P: usize = 5;
type W = [u64; P];

#[derive(Default)]
struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
}
type Map<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;

#[derive(Clone, Copy, PartialEq)]
enum Model {
    Nest,
    Stack,
}

/// Moves for one labelled piece: (slot index that supplies the multiplicity or None for the table,
/// consumed vector, pending vector). Stack state (s, f, e, t); nest state (c, t).
fn moves(model: Model, kind: usize) -> Vec<(Option<usize>, [i8; 4], [i8; 4])> {
    match model {
        Model::Stack => match kind {
            0 => vec![
                (None, [0, 0, 0, 0], [1, 0, 1, 0]),
                (Some(0), [1, 0, 0, 0], [1, 0, 1, 0]),
                (Some(1), [0, 1, 0, 0], [1, 0, 1, 0]),
                (Some(2), [0, 0, 1, 0], [0, 0, 1, 0]),
            ],
            1 => vec![
                (None, [0, 0, 0, 0], [0, 1, 0, 0]),
                (Some(0), [1, 0, 0, 0], [0, 1, 0, 0]),
                (Some(1), [0, 1, 0, 0], [0, 1, 0, 0]),
                (Some(2), [0, 0, 1, 0], [0, 0, 1, 0]),
            ],
            _ => vec![
                (None, [0, 0, 0, 0], [1, 0, 0, 1]),
                (Some(0), [1, 0, 0, 0], [1, 0, 0, 1]),
                (Some(3), [0, 0, 0, 1], [0, 0, 0, 1]),
            ],
        },
        Model::Nest => match kind {
            0 | 1 => vec![(None, [0, 0, 0, 0], [1, 0, 0, 0]), (Some(0), [1, 0, 0, 0], [1, 0, 0, 0])],
            _ => vec![(None, [0, 0, 0, 0], [0, 1, 0, 0]), (Some(1), [0, 1, 0, 0], [0, 1, 0, 0])],
        },
    }
}

struct Ctx {
    model: Model,
    width: usize,
    primes: [u64; P],
    mv: [Vec<(Option<usize>, [i8; 4], [i8; 4])>; 3],
}

fn wmul(a: &W, b: &W, pr: &[u64; P]) -> W {
    let mut r = [0; P];
    for i in 0..P {
        r[i] = mulmod(a[i], b[i], pr[i]);
    }
    r
}
fn wadd(a: &mut W, b: &W, pr: &[u64; P]) {
    for i in 0..P {
        let s = a[i] + b[i];
        a[i] = if s >= pr[i] { s - pr[i] } else { s };
    }
}
fn wsmall(k: u64, pr: &[u64; P]) -> W {
    let mut r = [0; P];
    for i in 0..P {
        r[i] = k % pr[i];
    }
    r
}

/// Place nD, nB, nT labelled pieces (L-cycles) of one size into slot state `st`.
fn place_group(ctx: &Ctx, st: &[u8], n: [u8; 3], l: u64) -> Vec<(Vec<u8>, W)> {
    let w = ctx.width;
    let mut cur: Map<(Vec<u8>, Vec<u8>), W> = Map::default();
    cur.insert((st.to_vec(), vec![0; w]), wsmall(1, &ctx.primes));
    for kind in 0..3 {
        for _ in 0..n[kind] {
            let mut nxt: Map<(Vec<u8>, Vec<u8>), W> = Map::default();
            for ((s, pend), val) in cur.iter() {
                for (key, cons, add) in &ctx.mv[kind] {
                    let mult = match key {
                        None => 1,
                        Some(i) => s[*i] as u64 * l,
                    };
                    if mult == 0 {
                        continue;
                    }
                    let mut s2 = s.clone();
                    let mut p2 = pend.clone();
                    for i in 0..w {
                        s2[i] = (s2[i] as i16 - cons[i] as i16) as u8;
                        p2[i] = (p2[i] as i16 + add[i] as i16) as u8;
                    }
                    let v = wmul(val, &wsmall(mult, &ctx.primes), &ctx.primes);
                    wadd(cur_entry(&mut nxt, (s2, p2)), &v, &ctx.primes);
                }
            }
            cur = nxt;
        }
    }
    let mut out: Map<Vec<u8>, W> = Map::default();
    for ((s, pend), val) in cur {
        let s2: Vec<u8> = s.iter().zip(&pend).map(|(a, b)| a + b).collect();
        wadd(cur_entry(&mut out, s2), &val, &ctx.primes);
    }
    out.into_iter().collect()
}

fn cur_entry<'a, K: std::hash::Hash + Eq>(m: &'a mut Map<K, W>, k: K) -> &'a mut W {
    m.entry(k).or_insert([0; P])
}

fn partitions(n: usize, maxp: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut out = vec![];
    for p in (1..=n.min(maxp)).rev() {
        for mut r in partitions(n - p, p) {
            r.insert(0, p);
            out.push(r);
        }
    }
    out
}

/// z_lambda mod p
fn zmod(lam: &[usize], p: u64) -> u64 {
    let mut r = 1u64;
    let mut i = 0;
    while i < lam.len() {
        let mut j = i;
        while j < lam.len() && lam[j] == lam[i] {
            j += 1;
        }
        let k = (j - i) as u64;
        r = mulmod(r, powmod(lam[i] as u64, k, p), p);
        for x in 1..=k {
            r = mulmod(r, x, p);
        }
        i = j;
    }
    r
}

struct Choice {
    w: W,
    atoms: Vec<(usize, [u8; 3])>, // (L, counts of D, B, T cycles of length L)
}

fn choices(m: usize, primes: &[u64; P]) -> Vec<Choice> {
    let mut out = vec![];
    for o in 0..=m {
        let c = m - o;
        for ld in partitions(c, c) {
            for lb in partitions(o, o) {
                for lt in partitions(o, o) {
                    let mut w = [0; P];
                    for i in 0..P {
                        let p = primes[i];
                        let z = mulmod(mulmod(zmod(&ld, p), zmod(&lb, p), p), zmod(&lt, p), p);
                        w[i] = powmod(z, p - 2, p);
                    }
                    let mut atoms: Vec<(usize, [u8; 3])> = vec![];
                    for (kind, lam) in [&ld, &lb, &lt].iter().enumerate() {
                        for &l in lam.iter() {
                            match atoms.iter_mut().find(|(x, _)| *x == l) {
                                Some((_, cnt)) => cnt[kind] += 1,
                                None => {
                                    let mut cnt = [0u8; 3];
                                    cnt[kind] = 1;
                                    atoms.push((l, cnt));
                                }
                            }
                        }
                    }
                    atoms.sort();
                    out.push(Choice { w, atoms });
                }
            }
        }
    }
    out
}

fn main() {
    disable_power_throttling();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let model = match args[0].as_str() {
        "nest" => Model::Nest,
        "stack" => Model::Stack,
        _ => panic!("model must be nest or stack"),
    };
    let n: usize = args[1].parse().unwrap();
    let width = if model == Model::Stack { 4 } else { 2 };
    let pv = primes_below(1u64 << 61, P);
    let mut primes = [0u64; P];
    primes.copy_from_slice(&pv);
    let mv = [moves(model, 0), moves(model, 1), moves(model, 2)];
    let mv = mv.map(|v| v.into_iter().map(|(k, c, a)| (k, c, a)).collect::<Vec<_>>());
    let ctx = Ctx { model, width, primes, mv };
    let _ = ctx.model;
    let ch: Vec<Vec<Choice>> = (0..=n).map(|m| if m == 0 { vec![] } else { choices(m, &primes) }).collect();
    let nthreads = threads();
    println!("# primes {}", primes.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(" "));

    // by_used[u]: map from state (n * width bytes: slot state per cycle length L = 1..n) to weight
    let mut by_used: Vec<Map<Vec<u8>, W>> = (0..=n).map(|_| Map::default()).collect();
    by_used[0].insert(vec![0u8; n * width], wsmall(1, &primes));
    let t0 = Instant::now();
    for used in 0..n {
        let ts = Instant::now();
        let level: Vec<(Vec<u8>, W)> = std::mem::take(&mut by_used[used]).into_iter().collect();
        if used > 0 {
            let mut tot = [0u64; P];
            for (_, v) in &level {
                wadd(&mut tot, v, &primes);
            }
            println!("{used} {}", tot.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "));
        }
        let nstates = level.len();
        // split work over threads; each returns its own maps for targets used+1..=n
        let chunk = (level.len() + nthreads - 1) / nthreads.max(1);
        let results: Vec<Vec<Map<Vec<u8>, W>>> = std::thread::scope(|sc| {
            let handles: Vec<_> = level
                .chunks(chunk.max(1))
                .map(|part| {
                    let ctx = &ctx;
                    let ch = &ch;
                    sc.spawn(move || {
                        let mut memo: Map<(Vec<u8>, [u8; 3], usize), Vec<(Vec<u8>, W)>> = Map::default();
                        let mut out: Vec<Map<Vec<u8>, W>> = (0..=n).map(|_| Map::default()).collect();
                        for (state, wt) in part {
                            for m in 1..=n - used {
                                for c in &ch[m] {
                                    let mut combos: Vec<(Vec<u8>, W)> = vec![(state.clone(), wmul(wt, &c.w, &ctx.primes))];
                                    for &(l, cnt) in &c.atoms {
                                        let seg = &state[(l - 1) * width..l * width];
                                        let key = (seg.to_vec(), cnt, l);
                                        if !memo.contains_key(&key) {
                                            let r = place_group(ctx, seg, cnt, l as u64);
                                            memo.insert(key.clone(), r);
                                        }
                                        let res = &memo[&key];
                                        let mut nxt = Vec::with_capacity(combos.len() * res.len());
                                        for (k0, w0) in &combos {
                                            for (ns, w1) in res {
                                                let mut k1 = k0.clone();
                                                k1[(l - 1) * width..l * width].copy_from_slice(ns);
                                                nxt.push((k1, wmul(w0, w1, &ctx.primes)));
                                            }
                                        }
                                        combos = nxt;
                                    }
                                    let tgt = &mut out[used + m];
                                    for (k1, w1) in combos {
                                        wadd(cur_entry(tgt, k1), &w1, &ctx.primes);
                                    }
                                }
                            }
                        }
                        out
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        for res in results {
            for (u, mp) in res.into_iter().enumerate() {
                if mp.is_empty() {
                    continue;
                }
                let tgt = &mut by_used[u];
                for (k, v) in mp {
                    wadd(cur_entry(tgt, k), &v, &primes);
                }
            }
        }
        let pending: usize = by_used.iter().map(|m| m.len()).sum();
        eprintln!(
            "level {used}: {nstates} states, {:.1}s (total {:.1}s), {pending} pending states",
            ts.elapsed().as_secs_f64(),
            t0.elapsed().as_secs_f64()
        );
    }
    let mut tot = [0u64; P];
    for v in by_used[n].values() {
        wadd(&mut tot, v, &primes);
    }
    println!("{n} {}", tot.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "));
    eprintln!("level {n} final, total {:.1}s", t0.elapsed().as_secs_f64());
}
