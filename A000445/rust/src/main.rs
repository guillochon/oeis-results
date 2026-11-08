// Consecutive k-th power residues: decide D(k, N):
//   is there a character chi: primes -> Z/k (completely additive "index mod k"),
//   with chi(2) in an allowed set, such that no r in [1, N-1] has chi(r) = chi(r+1) = 0 ?
// a(k) = max N with D(k, N) satisfiable (restricted to the residue classes of p that
// realize the order-k character; other p reduce to smaller k).
//
// Search: DFS over primes in increasing order. Assigning prime q determines chi on every
// n with largest prime factor q. Forward checking: a number whose largest prime Q is
// unassigned but whose cofactor is determined is "one-unassigned"; if a neighbor is a
// determined zero, values of chi(Q) making it zero are removed from dom[Q].
// Free-prime cutoff: once q > sqrt(N), all multiples m*q <= N have m < q determined; if
// chi({1..floor(N/q)}) misses some value v, set chi(q') = -v-ish for every prime q' >= q
// (each has fewer multiples), making all their multiples non-residues -> SAT.
// Symmetry: chi -> u*chi for units u of Z/k preserves the zero set; orbit-min pruning.

use std::env;
use std::io::Write as IoWrite;
use std::time::Instant;

struct Ctx {
    k: usize,
    n: usize,
    primes: Vec<u32>,
    pidx: Vec<u32>,    // prime -> index (only valid on primes)
    p1i: Vec<u32>,     // index of largest prime factor (u32::MAX for 1)
    e1: Vec<u8>,       // its exponent
    cof: Vec<u32>,     // n / P1^e1
    p2i: Vec<u32>,     // index of largest prime of cof (u32::MAX if cof == 1)
    by_p1_off: Vec<u32>,
    by_p1: Vec<u32>,
    by_p2_off: Vec<u32>,
    by_p2: Vec<u32>,
    chi: Vec<u8>,
    val: Vec<u8>,
    dom: Vec<u16>,
    trail: Vec<(u32, u16)>,
    units: Vec<usize>,
    nodes: Vec<u64>,
    total_nodes: u64,
    sat_depth: usize,
    max_depth: usize,
    t0: Instant,
    node_limit: u64,
    aborted: bool,
    reason: Vec<u32>,
    conf: (u8, u32), // (1, r) pair conflict; (2, qi) wipeout
    cert: Option<std::io::BufWriter<std::fs::File>>,
    two_mask: u16,
}

impl Ctx {
    fn new(k: usize, n: usize, two_mask: u16) -> Ctx {
        let mut lpf = vec![0u32; n + 1];
        let mut primes = Vec::new();
        for p in 2..=n {
            if lpf[p] == 0 {
                primes.push(p as u32);
                let mut m = p;
                while m <= n {
                    lpf[m] = p as u32;
                    m += p;
                }
            }
        }
        let mut pidx = vec![u32::MAX; n + 1];
        for (i, &p) in primes.iter().enumerate() {
            pidx[p as usize] = i as u32;
        }
        let mut p1i = vec![u32::MAX; n + 1];
        let mut e1 = vec![0u8; n + 1];
        let mut cof = vec![1u32; n + 1];
        let mut p2i = vec![u32::MAX; n + 1];
        for m in 2..=n {
            let p = lpf[m] as usize;
            let mut c = m;
            let mut e = 0u8;
            while c % p == 0 {
                c /= p;
                e += 1;
            }
            p1i[m] = pidx[p];
            e1[m] = e;
            cof[m] = c as u32;
            p2i[m] = if c == 1 { u32::MAX } else { pidx[lpf[c] as usize] };
        }
        let np = primes.len();
        let build = |key: &Vec<u32>| -> (Vec<u32>, Vec<u32>) {
            let mut cnt = vec![0u32; np + 1];
            for m in 2..=n {
                if key[m] != u32::MAX {
                    cnt[key[m] as usize + 1] += 1;
                }
            }
            for i in 0..np {
                cnt[i + 1] += cnt[i];
            }
            let mut pos = cnt.clone();
            let mut arr = vec![0u32; cnt[np] as usize];
            for m in 2..=n {
                if key[m] != u32::MAX {
                    let i = key[m] as usize;
                    arr[pos[i] as usize] = m as u32;
                    pos[i] += 1;
                }
            }
            (cnt, arr)
        };
        let (by_p1_off, by_p1) = build(&p1i);
        let (by_p2_off, by_p2) = build(&p2i);
        let full: u16 = ((1u32 << k) - 1) as u16;
        let mut dom = vec![full; np];
        dom[0] &= two_mask;
        let units: Vec<usize> = (1..k).filter(|&u| gcd(u, k) == 1).collect();
        let mut chi = vec![0u8; n + 2];
        chi[1] = 0;
        let mut c = Ctx {
            k,
            n,
            primes,
            pidx,
            p1i,
            e1,
            cof,
            p2i,
            by_p1_off,
            by_p1,
            by_p2_off,
            by_p2,
            chi,
            val: vec![0u8; np],
            dom,
            trail: Vec::new(),
            units,
            nodes: vec![0; np + 1],
            total_nodes: 0,
            sat_depth: 0,
            max_depth: 0,
            t0: Instant::now(),
            node_limit: u64::MAX,
            aborted: false,
            reason: vec![0u32; np * k],
            conf: (0, 0),
            cert: None,
            two_mask,
        };
        // initial: 1 is a determined zero; prime powers adjacent to 1 (i.e. 2) are one-unassigned.
        if n >= 2 {
            c.forbid_zero_of(2, 1);
        }
        c
    }

    #[inline]
    fn determined(&self, m: usize, cur: u32) -> bool {
        m == 1 || self.p1i[m] <= cur
    }

    // m is one-unassigned (cofactor determined); remove values of chi(P1(m)) that make chi(m)=0.
    // returns false on wipeout
    #[inline]
    fn forbid_zero_of(&mut self, m: usize, zn: usize) -> bool {
        let k = self.k;
        let qi = self.p1i[m] as usize;
        let s = self.chi[self.cof[m] as usize] as usize;
        let e = self.e1[m] as usize;
        let mut bad: u16 = 0;
        for w in 0..k {
            if (s + e * w) % k == 0 {
                bad |= 1 << w;
            }
        }
        let d = self.dom[qi];
        if d & bad != 0 {
            let r = m.min(zn) as u32;
            for w in 0..k {
                if (d & bad) >> w & 1 == 1 {
                    self.reason[qi * k + w] = r;
                }
            }
            self.trail.push((qi as u32, d));
            self.dom[qi] = d & !bad;
        }
        if self.dom[qi] == 0 {
            self.conf = (2, qi as u32);
            return false;
        }
        true
    }

    fn undo(&mut self, mark: usize) {
        while self.trail.len() > mark {
            let (i, d) = self.trail.pop().unwrap();
            self.dom[i as usize] = d;
        }
    }

    // assign prime index idx = v; returns false on conflict
    fn assign(&mut self, idx: usize, v: usize) -> bool {
        let k = self.k;
        let n = self.n;
        let cur = idx as u32;
        self.val[idx] = v as u8;
        let (a, b) = (self.by_p1_off[idx] as usize, self.by_p1_off[idx + 1] as usize);
        for t in a..b {
            let m = self.by_p1[t] as usize;
            let s = self.chi[self.cof[m] as usize] as usize;
            self.chi[m] = ((s + self.e1[m] as usize * v) % k) as u8;
        }
        for t in a..b {
            let m = self.by_p1[t] as usize;
            if self.chi[m] != 0 {
                continue;
            }
            for nb in [m - 1, m + 1] {
                if nb < 1 || nb > n {
                    continue;
                }
                if self.determined(nb, cur) {
                    if self.chi[nb] == 0 {
                        self.conf = (1, m.min(nb) as u32);
                        return false;
                    }
                } else if self.p2i[nb] == u32::MAX || self.p2i[nb] <= cur {
                    if !self.forbid_zero_of(nb, m) {
                        return false;
                    }
                }
            }
        }
        let (a, b) = (self.by_p2_off[idx] as usize, self.by_p2_off[idx + 1] as usize);
        for t in a..b {
            let m = self.by_p2[t] as usize;
            for nb in [m - 1, m + 1] {
                if nb < 1 || nb > n {
                    continue;
                }
                if self.determined(nb, cur) && self.chi[nb] == 0 {
                    if !self.forbid_zero_of(m, nb) {
                        return false;
                    }
                    break;
                }
            }
        }
        true
    }

    // all primes from idx on can be made "all multiples nonzero"?
    fn free_from(&self, idx: usize) -> bool {
        let q = self.primes[idx] as usize;
        let mq = self.n / q;
        if mq < self.k {
            return true;
        }
        if q * q <= self.n {
            return false;
        }
        let mut seen: u32 = 0;
        for m in 1..=mq {
            seen |= 1 << self.chi[m];
        }
        seen != (1u32 << self.k) - 1
    }

    fn dfs(&mut self, idx: usize, group: &[usize]) -> bool {
        self.total_nodes += 1;
        self.nodes[idx] += 1;
        if idx > self.max_depth {
            self.max_depth = idx;
        }
        if self.total_nodes >= self.node_limit {
            self.aborted = true;
            return false;
        }
        if self.total_nodes % 50_000_000 == 0 {
            eprintln!(
                "  [{:.0}s] nodes={} max_depth={} (p={})",
                self.t0.elapsed().as_secs_f64(),
                self.total_nodes,
                self.max_depth,
                self.primes[self.max_depth.min(self.primes.len() - 1)]
            );
        }
        if idx >= self.primes.len() || self.free_from(idx) {
            self.sat_depth = idx;
            return true;
        }
        let k = self.k;
        let d = self.dom[idx];
        let base: u16 = if idx == 0 { self.two_mask } else { ((1u32 << k) - 1) as u16 };
        // orbit representatives of allowed values; value order: odd-ish first, 0 last
        let mut reps: Vec<usize> = (0..k).filter(|&v| base >> v & 1 == 1).collect();
        reps.retain(|&v| group.iter().all(|&u| (u * v) % k >= v));
        reps.sort_by_key(|&v| (gcd(v, k), v));
        if let Some(w) = self.cert.as_mut() {
            writeln!(w, "B {} {}", self.primes[idx], reps.len()).unwrap();
        }
        for v in reps {
            if d >> v & 1 == 0 {
                let r = self.reason[idx * k + v];
                if let Some(w) = self.cert.as_mut() {
                    writeln!(w, "V {}
X {}", v, r).unwrap();
                }
                continue;
            }
            let mark = self.trail.len();
            if let Some(w) = self.cert.as_mut() {
                writeln!(w, "V {}", v).unwrap();
            }
            let ok = self.assign(idx, v);
            if !ok && self.cert.is_some() {
                let (kind, x) = self.conf;
                let mut out = String::new();
                if kind == 1 {
                    out.push_str(&format!("X {}
", x));
                } else {
                    let qi = x as usize;
                    out.push_str(&format!("B {} {}
", self.primes[qi], k));
                    for w2 in 0..k {
                        out.push_str(&format!("V {}
X {}
", w2, self.reason[qi * k + w2]));
                    }
                }
                self.cert.as_mut().unwrap().write_all(out.as_bytes()).unwrap();
            }
            if ok {
                let g2: Vec<usize> = group.iter().cloned().filter(|&u| (u * v) % k == v).collect();
                if self.dfs(idx + 1, &g2) {
                    return true;
                }
                if self.aborted {
                    return false;
                }
            }
            self.undo(mark);
        }
        false
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn solve(k: usize, n: usize, two_mask: u16, limit: u64, verbose: bool, cert: Option<String>) -> (Option<bool>, Ctx) {
    let mut c = Ctx::new(k, n, two_mask);
    if let Some(path) = cert {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        writeln!(w, "K {} N {} TWOMASK {:b}", k, n, two_mask).unwrap();
        c.cert = Some(w);
    }
    c.node_limit = limit;
    let units = c.units.clone();
    let mut g = vec![1usize];
    g.extend(units.iter().filter(|&&u| u != 1));
    let r = c.dfs(0, &g);
    if verbose {
        eprintln!(
            "k={} N={} -> {} nodes={} maxdepth={} ({:.2}s)",
            k,
            n,
            if c.aborted { "ABORT" } else if r { "SAT" } else { "UNSAT" },
            c.total_nodes,
            c.max_depth,
            c.t0.elapsed().as_secs_f64()
        );
    }
    (if c.aborted { None } else { Some(r) }, c)
}

fn main() {
    let h = std::thread::Builder::new().stack_size(1 << 30).spawn(real_main).unwrap();
    h.join().unwrap();
}

fn real_main() {
    let args: Vec<String> = env::args().collect();
    let mode = args[1].as_str();
    let k: usize = args[2].parse().unwrap();
    let two_mask: u16 = if k == 8 { 0b0101_0101 } else { ((1u32 << k) - 1) as u16 };
    let two_mask = match env::var("TWOMASK") {
        Ok(s) => u16::from_str_radix(&s, 2).unwrap(),
        Err(_) => two_mask,
    };
    let limit: u64 = env::var("NODELIMIT").ok().map(|s| s.parse().unwrap()).unwrap_or(u64::MAX);
    match mode {
        "decide" => {
            let n: usize = args[3].parse().unwrap();
            let (r, c) = solve(k, n, two_mask, limit, true, env::var("CERT").ok());
            if r == Some(true) {
                // print witness chi on primes up to sat_depth
                let w: Vec<String> = (0..c.sat_depth)
                    .map(|i| format!("{}:{}", c.primes[i], c.val[i]))
                    .collect();
                println!("SAT {} {}", n, c.sat_depth);
                if let Some(path) = args.get(4) {
                    std::fs::write(path, w.join("\n")).unwrap();
                }
            } else if r == Some(false) {
                println!("UNSAT {}", n);
                let nz: Vec<String> = c.nodes.iter().enumerate().filter(|(_, &x)| x > 0)
                    .map(|(i, x)| format!("{}:{}", c.primes[i.min(c.primes.len()-1)], x)).collect();
                eprintln!("nodes by prime: {}", nz.join(" "));
            } else {
                println!("ABORT {}", n);
            }
        }
        "max" => {
            // binary search for max N with D(N) SAT, in [lo, hi]
            let mut lo: usize = args[3].parse().unwrap(); // known SAT
            let mut hi: usize = args[4].parse().unwrap(); // known UNSAT
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                let (r, _) = solve(k, mid, two_mask, limit, true, None);
                match r {
                    Some(true) => lo = mid,
                    Some(false) => hi = mid,
                    None => {
                        println!("ABORT at {}", mid);
                        return;
                    }
                }
            }
            println!("k={} a(k)={}", k, lo);
        }
        _ => panic!("mode"),
    }
}
