//! Large-grid engine for A350148: the same transfer matrix as `main.rs` (it calls its `step`), with
//!   * 64-bit keys: strand ends stored as brackets (2 bits per frontier point, in order along the frontier),
//!     which works because strand pairings below a frontier never cross;
//!   * counts as 128-bit integers stored in two u64 words (24-byte records);
//!   * the next layer split into S hash shards built G groups at a time, so only 1/G of it is in memory;
//!     the current layer is kept in memory or spilled to disk when larger than the budget;
//!   * T worker threads.
//!
//! usage: big FILE [--threads T] [--mem-gb M] [--tmp DIR]

#[path = "../main.rs"]
#[allow(dead_code)]
mod core;

use core::{initial, read_problem, step, Problem, St, DEG2, MAXP};
use std::collections::HashMap;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

#[derive(Default)]
struct Mix64(u64);
impl Hasher for Mix64 {
    fn finish(&self) -> u64 {
        mix(self.0)
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!()
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = v;
    }
}
fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// A count, as two u64 words so a record is 24 bytes.
#[derive(Clone, Copy, Default)]
struct V([u64; 2]);
impl V {
    fn add(&mut self, o: V) {
        let (lo, c) = self.0[0].overflowing_add(o.0[0]);
        self.0[0] = lo;
        self.0[1] = self.0[1].wrapping_add(o.0[1]).wrapping_add(c as u64);
    }
    fn get(&self) -> u128 {
        (self.0[1] as u128) << 64 | self.0[0] as u128
    }
}
type Map = HashMap<u64, V, BuildHasherDefault<Mix64>>;

/// Slots in frontier order before cell column x: p[0..x], k, p[x..=n].
fn order(x: usize, n: usize) -> impl Iterator<Item = usize> {
    // 0..x -> p[i]; x -> k (index MAXP as a marker); then p[x..=n]
    (0..x).chain(std::iter::once(MAXP)).chain(x..=n)
}

fn pack(s: &St, x: usize) -> u64 {
    let n = s.n as usize;
    let mut seen = [false; 16];
    let mut key: u64 = 0;
    for i in order(x, n) {
        let c = if i == MAXP { s.k } else { s.p[i] };
        let code = if c == 0 {
            0
        } else if c == DEG2 {
            3
        } else if !seen[c as usize] {
            seen[c as usize] = true;
            1
        } else {
            2
        };
        key = key << 2 | code;
    }
    key = key << n | (s.h as u64 & ((1u64 << n) - 1));
    key << 1 | s.v as u64
}

fn unpack(mut key: u64, n: usize, x: usize) -> St {
    let v = key & 1 == 1;
    key >>= 1;
    let h = (key & ((1u64 << n) - 1)) as u32;
    key >>= n;
    let mut slots = [0usize; MAXP + 2];
    for (j, i) in order(x, n).enumerate() {
        slots[j] = i;
    }
    let ns = n + 2;
    let mut codes = [0u8; MAXP + 2];
    for j in (0..ns).rev() {
        codes[j] = (key & 3) as u8;
        key >>= 2;
    }
    let mut s = St { p: [0u8; MAXP], k: 0, h, v, n: n as u8 };
    let mut stack = [0u8; 16];
    let mut sp = 0usize;
    let mut next = 1u8;
    for (j, &i) in slots[..ns].iter().enumerate() {
        let c = match codes[j] {
            0 => 0,
            3 => DEG2,
            1 => {
                stack[sp] = next;
                sp += 1;
                next += 1;
                next - 1
            }
            _ => {
                sp -= 1;
                stack[sp]
            }
        };
        if i == MAXP {
            s.k = c;
        } else {
            s.p[i] = c;
        }
    }
    s
}

const S: usize = 256; // shards
const BUF: usize = 4096; // records per thread-local shard buffer

enum Chunk {
    Mem(Vec<(u64, V)>),
    Disk(PathBuf, usize),
}

fn write_chunk(path: &PathBuf, recs: &[(u64, V)]) {
    let mut w = BufWriter::with_capacity(1 << 20, File::create(path).expect("create chunk"));
    for (k, v) in recs {
        w.write_all(&k.to_le_bytes()).unwrap();
        w.write_all(&v.0[0].to_le_bytes()).unwrap();
        w.write_all(&v.0[1].to_le_bytes()).unwrap();
    }
}

fn read_chunk(path: &PathBuf, len: usize) -> Vec<(u64, V)> {
    let mut r = BufReader::with_capacity(1 << 20, File::open(path).expect("open chunk"));
    let mut out = Vec::with_capacity(len);
    let mut b = [0u8; 24];
    for _ in 0..len {
        r.read_exact(&mut b).unwrap();
        let k = u64::from_le_bytes(b[0..8].try_into().unwrap());
        let lo = u64::from_le_bytes(b[8..16].try_into().unwrap());
        let hi = u64::from_le_bytes(b[16..24].try_into().unwrap());
        out.push((k, V([lo, hi])));
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pr: Problem = read_problem(&args[1]);
    let n = pr.n;
    assert!(2 * (n + 2) + n + 1 <= 64, "key does not fit in 64 bits");
    let mut threads = 6usize;
    let mut mem_gb = 6.0f64;
    let mut tmp = std::env::temp_dir().join(format!("a350148_{}", std::process::id()));
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--threads" => threads = args[i + 1].parse().unwrap(),
            "--mem-gb" => mem_gb = args[i + 1].parse().unwrap(),
            "--tmp" => tmp = PathBuf::from(&args[i + 1]),
            a => panic!("unknown option {}", a),
        }
        i += 2;
    }
    std::fs::create_dir_all(&tmp).unwrap();
    let budget = (mem_gb * 1e9) as usize;
    const MAP_BYTES: usize = 36; // per entry in a hash map, including load-factor slack
    const REC: usize = 24;

    let t0 = Instant::now();
    let s0 = initial(n);
    let mut cur: Vec<Chunk> = vec![Chunk::Mem(vec![(pack(&s0, 0), V([1, 0]))])];
    let mut cur_len = 1usize;
    let mut peak = 0usize;
    for idx in 0..n * n {
        let x = idx % n;
        let x2 = (idx + 1) % n;
        // estimate the next layer as 1.5x this one; choose the group count to fit the budget
        let est = cur_len * 3 / 2 + 1;
        let cur_mem: usize = cur.iter().map(|c| if let Chunk::Mem(v) = c { v.len() * REC } else { 0 }).sum();
        let avail = budget.saturating_sub(cur_mem).max(budget / 4);
        let mut g = 1usize;
        while g < S && est * MAP_BYTES / g > avail {
            g *= 2;
        }
        let mut next: Vec<Option<Chunk>> = (0..S).map(|_| None).collect();
        let mut next_len = 0usize;
        for grp in 0..g {
            // pre-size each shard map for about one current layer's worth of states, so it rarely rehashes
            let cap = cur_len / S + 16;
            let maps: Vec<Mutex<Map>> =
                (0..S / g).map(|_| Mutex::new(Map::with_capacity_and_hasher(cap, Default::default()))).collect();
            let work = AtomicUsize::new(0);
            std::thread::scope(|sc| {
                for _ in 0..threads {
                    sc.spawn(|| {
                        let mut bufs: Vec<Vec<(u64, V)>> = (0..S / g).map(|_| Vec::with_capacity(BUF)).collect();
                        let mut succ = Vec::new();
                        let flush = |b: &mut Vec<(u64, V)>, m: &Mutex<Map>| {
                            let mut m = m.lock().unwrap();
                            for (k, v) in b.drain(..) {
                                m.entry(k).or_default().add(v);
                            }
                        };
                        loop {
                            let ci = work.fetch_add(1, Ordering::Relaxed);
                            if ci >= cur.len() {
                                break;
                            }
                            let loaded;
                            let recs: &[(u64, V)] = match &cur[ci] {
                                Chunk::Mem(v) => v,
                                Chunk::Disk(p, len) => {
                                    loaded = read_chunk(p, *len);
                                    &loaded
                                }
                            };
                            for &(key, val) in recs {
                                let s = unpack(key, n, x);
                                succ.clear();
                                step(&pr, idx, &s, &mut succ);
                                for (_, ns) in succ.drain(..) {
                                    let k2 = pack(&ns, x2);
                                    let sh = (mix(k2) >> 32) as usize % S;
                                    if sh % g != grp {
                                        continue;
                                    }
                                    let b = &mut bufs[sh / g];
                                    b.push((k2, val));
                                    if b.len() >= BUF {
                                        flush(b, &maps[sh / g]);
                                    }
                                }
                            }
                        }
                        for (j, b) in bufs.iter_mut().enumerate() {
                            if !b.is_empty() {
                                flush(b, &maps[j]);
                            }
                        }
                    });
                }
            });
            // turn this group's maps into chunks: in memory while it fits, else on disk
            for (j, m) in maps.into_iter().enumerate() {
                let m = m.into_inner().unwrap();
                let sh = j * g + grp;
                next_len += m.len();
                let recs: Vec<(u64, V)> = m.into_iter().collect();
                let spill = g > 1 || (cur_len + next_len) * REC > budget;
                next[sh] = Some(if spill && !recs.is_empty() {
                    let p = tmp.join(format!("l{}_s{}.bin", (idx + 1) % 2, sh));
                    write_chunk(&p, &recs);
                    Chunk::Disk(p, recs.len())
                } else {
                    Chunk::Mem(recs)
                });
            }
        }
        // the old layer's files are overwritten two layers later; delete them now to free disk
        for c in &cur {
            if let Chunk::Disk(p, _) = c {
                let _ = std::fs::remove_file(p);
            }
        }
        cur = next.into_iter().flatten().collect();
        cur_len = next_len;
        peak = peak.max(cur_len);
        if (idx + 1) % n == 0 {
            eprintln!("row {:2}/{}  states {:>11}  groups {:>3}  {:>8.1}s", (idx + 1) / n, n, cur_len, g,
                      t0.elapsed().as_secs_f64());
        }
    }
    let mut total = V::default();
    for c in &cur {
        let recs = match c {
            Chunk::Mem(v) => v.clone(),
            Chunk::Disk(p, len) => read_chunk(p, *len),
        };
        for (_, v) in recs {
            total.add(v);
        }
    }
    for c in &cur {
        if let Chunk::Disk(p, _) = c {
            let _ = std::fs::remove_file(p);
        }
    }
    let _ = std::fs::remove_dir(&tmp);
    // same result line as `rust count` (count, then the peak and time)
    println!("peak states {} {} n={} time {:.1}s", peak, total.get(), n, t0.elapsed().as_secs_f64());
}
