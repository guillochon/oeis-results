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
    n: usize,      // largest number considered = nwin + l - 2
    l: usize,      // window length (2 = pairs, 3 = triples)
    nwin: usize,   // windows start at r in [1, nwin-1]
    last_beat: f64,
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
    hint: Vec<u8>,   // preferred value per prime index (255 = none); SAT search only
    rng: u64,        // 0 = deterministic value order; else xorshift state for shuffling
    deadline: f64,   // abort when elapsed seconds exceed this
    cbj: bool,       // conflict-directed backjumping (off when writing a certificate)
    cs: Vec<u32>,    // conflict set (prime indices) returned by a failed dfs call
    gshrink: Vec<u32>, // levels at which the symmetry group shrank (symmetry pruning depends on them)
}

impl Ctx {
    fn new(k: usize, nwin: usize, l: usize, two_mask: u16) -> Ctx {
        let n = nwin + l - 2;
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
        // by_p2 only needs numbers m that, at the moment they become one-unassigned (when their
        // second-largest prime p2i[m] is assigned), sit in some window whose other members are all
        // already determined (p1i <= p2i[m]). This is static, so filter once; propagation is unchanged.
        let mut key2 = p2i.clone();
        for m in 2..=n {
            let q = p2i[m];
            if q == u32::MAX {
                continue;
            }
            let lo = if m + 1 > l { m + 1 - l } else { 1 };
            let hi = m.min(nwin - 1);
            let mut keep = false;
            for st in lo..=hi {
                if (st..st + l).all(|x| x == m || x == 1 || (x <= n && p1i[x] <= q)) {
                    keep = true;
                    break;
                }
            }
            if !keep {
                key2[m] = u32::MAX;
            }
        }
        let (by_p2_off, by_p2) = build(&key2);
        drop(key2);
        let full: u16 = ((1u32 << k) - 1) as u16;
        let mut dom = vec![full; np];
        dom[0] &= two_mask;
        let units: Vec<usize> = (1..k).filter(|&u| gcd(u, k) == 1).collect();
        let mut chi = vec![0u8; n + 2];
        chi[1] = 0;
        let mut c = Ctx {
            k,
            n,
            l,
            nwin,
            last_beat: 0.0,
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
            hint: vec![255u8; np],
            rng: 0,
            deadline: f64::INFINITY,
            cbj: env::var("NOCBJ").is_err(),
            cs: Vec::new(),
            gshrink: Vec::new(),
        };
        // initial: 1 is a determined zero; prime powers adjacent to 1 (i.e. 2) are one-unassigned.
        if l == 2 && nwin >= 2 {
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
    fn forbid_zero_of(&mut self, m: usize, r: usize) -> bool {
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
            let r = r as u32;
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

    // prime indices involved in window st (all prime factors of its members)
    fn culprits(&self, st: usize, out: &mut Vec<u32>) {
        for x in st..st + self.l {
            let mut y = x;
            while y > 1 {
                out.push(self.p1i[y]);
                y = self.cof[y] as usize;
            }
        }
    }

    fn undo(&mut self, mark: usize) {
        while self.trail.len() > mark {
            let (i, d) = self.trail.pop().unwrap();
            self.dom[i as usize] = d;
        }
    }

    #[inline]
    fn det_zero(&self, x: usize, cur: u32) -> bool {
        self.determined(x, cur) && self.chi[x] == 0
    }

    // assign prime index idx = v; returns false on conflict.
    // A window is W(s) = {s, ..., s+l-1} with 1 <= s <= nwin-1; it is a conflict when all members are zero.
    fn assign(&mut self, idx: usize, v: usize) -> bool {
        let k = self.k;
        let l = self.l;
        let cur = idx as u32;
        self.val[idx] = v as u8;
        let (a, b) = (self.by_p1_off[idx] as usize, self.by_p1_off[idx + 1] as usize);
        for t in a..b {
            let m = self.by_p1[t] as usize;
            let s0 = self.chi[self.cof[m] as usize] as usize;
            self.chi[m] = ((s0 + self.e1[m] as usize * v) % k) as u8;
        }
        for t in a..b {
            let m = self.by_p1[t] as usize;
            if self.chi[m] != 0 {
                continue;
            }
            let lo = if m + 1 > l { m + 1 - l } else { 1 };
            let hi = m.min(self.nwin - 1);
            for st in lo..=hi {
                // members other than m: find the not-determined-zero ones
                let mut open = usize::MAX;
                let mut nopen = 0;
                let mut dead = false;
                for x in st..st + l {
                    if x == m || self.det_zero(x, cur) {
                        continue;
                    }
                    if self.determined(x, cur) {
                        dead = true;
                        break;
                    }
                    nopen += 1;
                    open = x;
                    if nopen > 1 {
                        break;
                    }
                }
                if dead || nopen > 1 {
                    continue;
                }
                if nopen == 0 {
                    self.conf = (1, st as u32);
                    return false;
                }
                if self.p2i[open] == u32::MAX || self.p2i[open] <= cur {
                    if !self.forbid_zero_of(open, st) {
                        return false;
                    }
                }
            }
        }
        let (a, b) = (self.by_p2_off[idx] as usize, self.by_p2_off[idx + 1] as usize);
        for t in a..b {
            let m = self.by_p2[t] as usize; // undetermined, now one-unassigned
            let lo = if m + 1 > l { m + 1 - l } else { 1 };
            let hi = m.min(self.nwin - 1);
            for st in lo..=hi {
                if (st..st + l).all(|x| x == m || self.det_zero(x, cur)) {
                    if !self.forbid_zero_of(m, st) {
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
        if self.total_nodes % 1024 == 0 && self.t0.elapsed().as_secs_f64() > self.deadline {
            self.aborted = true;
            return false;
        }
        if self.total_nodes % 1024 == 0 && self.t0.elapsed().as_secs_f64() - self.last_beat >= 30.0 {
            self.last_beat = self.t0.elapsed().as_secs_f64();
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
        let nbase = reps.len();
        reps.retain(|&v| group.iter().all(|&u| (u * v) % k >= v));
        let cbj = self.cbj && self.cert.is_none();
        let mut acc: Vec<u32> = Vec::new();
        if cbj && reps.len() < nbase {
            acc.extend(self.gshrink.iter().cloned());
        }
        reps.sort_by_key(|&v| (gcd(v, k), v));
        if self.cert.is_none() {
            if self.rng != 0 {
                // shuffle within each gcd class (keeps units before non-units, 0 last)
                let mut i = 0;
                while i < reps.len() {
                    let g = gcd(reps[i], k);
                    let mut j = i;
                    while j < reps.len() && gcd(reps[j], k) == g {
                        j += 1;
                    }
                    for t in (i + 1..j).rev() {
                        self.rng ^= self.rng << 13;
                        self.rng ^= self.rng >> 7;
                        self.rng ^= self.rng << 17;
                        let r = i + (self.rng % ((t - i + 1) as u64)) as usize;
                        reps.swap(t, r);
                    }
                    i = j;
                }
            }
            let h = self.hint[idx] as usize;
            if h < k && base >> h & 1 == 1 {
                // a hinted value need not be an orbit representative; it is still a valid
                // choice for SAT search (no UNSAT claims are made in this mode)
                reps.retain(|&v| v != h);
                reps.insert(0, h);
            }
        }
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
                if cbj {
                    self.culprits(r as usize, &mut acc);
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
            if !ok && cbj {
                let (kind, x) = self.conf;
                if kind == 1 {
                    self.culprits(x as usize, &mut acc);
                } else {
                    let qi = x as usize;
                    for w2 in 0..k {
                        let r = self.reason[qi * k + w2] as usize;
                        self.culprits(r, &mut acc);
                    }
                }
            }
            if ok {
                let g2: Vec<usize> = group.iter().cloned().filter(|&u| (u * v) % k == v).collect();
                let shrank = g2.len() < group.len();
                if shrank {
                    self.gshrink.push(idx as u32);
                }
                let res = self.dfs(idx + 1, &g2);
                if shrank {
                    self.gshrink.pop();
                }
                if res {
                    return true;
                }
                if self.aborted {
                    return false;
                }
                if cbj {
                    let cs = std::mem::take(&mut self.cs);
                    if !cs.contains(&(idx as u32)) {
                        // this level is not responsible: jump over it
                        self.undo(mark);
                        self.cs = cs;
                        return false;
                    }
                    acc.extend(cs);
                }
            }
            self.undo(mark);
        }
        if cbj {
            acc.retain(|&x| (x as usize) < idx);
            acc.sort_unstable();
            acc.dedup();
            self.cs = acc;
        }
        false
    }
}

fn luby(i: u64) -> u64 {
    // i >= 1
    let mut k = 1u32;
    while (1u64 << k) - 1 < i {
        k += 1;
    }
    if (1u64 << k) - 1 == i {
        1u64 << (k - 1)
    } else {
        luby(i - ((1u64 << (k - 1)) - 1))
    }
}

/// SAT-seeking with restarts: first pass follows the hint (previous witness) deterministically,
/// later passes shuffle value order with Luby-scheduled node limits. Returns witness or None on timeout.
fn seek(k: usize, n: usize, two_mask: u16, secs: f64, unit: u64, hint: Option<String>) -> Option<Vec<String>> {
    let l: usize = env::var("WINDOW").ok().map(|s| s.parse().unwrap()).unwrap_or(2);
    let mut c = Ctx::new(k, n, l, two_mask);
    eprintln!("  built tables for k={} N={} L={} in {:.1}s", k, n, l, c.t0.elapsed().as_secs_f64());
    if let Some(path) = hint {
        let mut cnt = 0;
        for line in std::fs::read_to_string(path).unwrap().lines() {
            if let Some((p, v)) = line.trim().split_once(':') {
                let p: usize = p.parse().unwrap();
                if p <= c.n && c.pidx[p] != u32::MAX {
                    c.hint[c.pidx[p] as usize] = (v.parse::<usize>().unwrap() % k) as u8;
                    cnt += 1;
                }
            }
        }
        eprintln!("  hint: {} primes", cnt);
    }
    c.deadline = secs;
    let mut g = vec![1usize];
    g.extend(c.units.iter().filter(|&&u| u != 1));
    let mark0 = c.trail.len();
    let mut i: u64 = 1;
    loop {
        let lim = luby(i) * unit;
        c.node_limit = c.total_nodes + lim;
        c.aborted = false;
        c.max_depth = 0;
        c.rng = if i == 1 { 0 } else { 0x9E3779B97F4A7C15u64.wrapping_mul(i) | 1 };
        let r = c.dfs(0, &g);
        let el = c.t0.elapsed().as_secs_f64();
        if r {
            eprintln!("k={} L={} N={} -> SAT on restart {} nodes={} ({:.2}s)", k, l, n, i, c.total_nodes, el);
            return Some((0..c.sat_depth).map(|j| format!("{}:{}", c.primes[j], c.val[j])).collect());
        }
        c.undo(mark0);
        if !c.aborted {
            // exhausted the whole tree without hint/shuffle-induced incompleteness only on pass 1
            if i == 1 {
                eprintln!("k={} L={} N={} -> UNSAT (complete pass) nodes={} ({:.2}s)", k, l, n, c.total_nodes, el);
                println!("UNSAT {}", n);
                std::process::exit(0);
            }
        }
        eprintln!("  restart {} limit={} maxdepth={} (p={}) total={} {:.0}s", i, lim, c.max_depth,
            c.primes[c.max_depth.min(c.primes.len() - 1)], c.total_nodes, el);
        if el > secs {
            eprintln!("k={} L={} N={} -> TIMEOUT nodes={} ({:.2}s)", k, l, n, c.total_nodes, el);
            return None;
        }
        i += 1;
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn solve(k: usize, n: usize, two_mask: u16, limit: u64, verbose: bool, cert: Option<String>) -> (Option<bool>, Ctx) {
    let l: usize = env::var("WINDOW").ok().map(|s| s.parse().unwrap()).unwrap_or(2);
    let mut c = Ctx::new(k, n, l, two_mask);
    if verbose {
        eprintln!("  built tables for k={} N={} L={} in {:.1}s", k, n, l, c.t0.elapsed().as_secs_f64());
    }
    if let Some(path) = cert {
        let mut w = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        writeln!(w, "K {} N {} TWOMASK {:b} L {}", k, n, two_mask, l).unwrap();
        c.cert = Some(w);
    }
    c.node_limit = limit;
    if let Ok(d) = env::var("DEADLINE") {
        c.deadline = d.parse().unwrap();
    }
    let units = c.units.clone();
    let mut g = vec![1usize];
    g.extend(units.iter().filter(|&&u| u != 1));
    let r = c.dfs(0, &g);
    if verbose {
        eprintln!(
            "k={} L={} N={} -> {} nodes={} maxdepth={} ({:.2}s)",
            k,
            c.l,
            n,
            if c.aborted { "ABORT" } else if r { "SAT" } else { "UNSAT" },
            c.total_nodes,
            c.max_depth,
            c.t0.elapsed().as_secs_f64()
        );
    }
    (if c.aborted { None } else { Some(r) }, c)
}

/// Windows 11 runs background console jobs as EcoQoS work on the E-cores; opt out.
#[cfg(windows)]
fn disable_power_throttling() {
    #[repr(C)]
    struct St {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
    }
    let st = St { version: 1, control_mask: 1, state_mask: 0 };
    let ok = unsafe {
        SetProcessInformation(GetCurrentProcess(), 4, &st as *const _ as *const std::ffi::c_void, std::mem::size_of::<St>() as u32)
    };
    if ok == 0 {
        eprintln!("note: could not disable Windows power throttling (EcoQoS)");
    }
}

#[cfg(not(windows))]
fn disable_power_throttling() {}

fn main() {
    disable_power_throttling();
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
        "seek" => {
            // seek k N secs witness_out [hint]
            let n: usize = args[3].parse().unwrap();
            let secs: f64 = args[4].parse().unwrap();
            let unit: u64 = env::var("UNIT").ok().map(|s| s.parse().unwrap()).unwrap_or(20000);
            match seek(k, n, two_mask, secs, unit, args.get(6).cloned()) {
                Some(w) => {
                    println!("SAT {} {}", n, w.len());
                    std::fs::write(&args[5], w.join("
")).unwrap();
                }
                None => println!("NOSAT {}", n),
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
