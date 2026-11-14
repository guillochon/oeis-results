//! Self-dual semigroups of order n by Burnside over anti-automorphisms.
//!
//!   SD(n) = (1/n!) sum_{pi in S_n} A(pi),   A(pi) = #{associative tables T on [n] : pi(xy) = pi(y) pi(x)}.
//!
//! For one representative pi per cycle type, A(pi)/|C(pi)| = sum over C(pi)-orbits of 1/|stabilizer|.
//! The orbits are enumerated by a backtracking search over the table with associativity propagation and
//! a lex-leader test against every element of C(pi) (kept as an explicit list that shrinks down the tree).
//! A leaf is the lex-least table of its orbit; the group elements still alive there are its stabilizer.
//!
//! Usage:
//!   selfdual sd <n> [threads] [first_class] [last_class]   self-dual count, all cycle types (or a range)
//!   selfdual iso <n> [threads]                             all semigroups up to isomorphism (no pi), a check
//! Output lines "class <i> <cycle type> k=<idempotents> nil3=<0/1> stab=<s> count=<c>"; the count of
//! self-dual semigroups with the given invariants is sum count/stab.

use std::collections::HashMap;
use std::time::Instant;

const MAXN: usize = 10;
const UNK: u8 = 255;
const NC: usize = MAXN * MAXN + MAXN; // table cells, then one depth cell per element
const INF: u8 = 15; // depth of the zero
const KIND_GROUP: u8 = 14; // depth-cell value of a unique idempotent that is not a zero

#[derive(Clone)]
struct Ctx {
    n: usize,
    anti: bool,
    pi: [u8; MAXN],
    order: Vec<u8>,          // cells in comparison / branching order: diagonal first, then row-major
    g: Vec<[u8; MAXN]>,      // the acting group
    ginv: Vec<[u8; MAXN]>,
    threads: usize,
    tid: usize,
    split: usize,
    no3nil: bool,
}

struct St {
    val: [u8; NC],
    dom: [u16; NC],
    trail: Vec<(u8, u16, u8)>,
    queue: Vec<u8>,
}

impl St {
    #[inline]
    fn set_dom(&mut self, c: usize, d: u16) {
        self.trail.push((c as u8, self.dom[c], self.val[c]));
        self.dom[c] = d;
    }
    #[inline]
    fn assign(&mut self, c: usize, v: u8) -> bool {
        if self.val[c] != UNK {
            return self.val[c] == v;
        }
        if self.dom[c] & (1 << v) == 0 {
            return false;
        }
        self.trail.push((c as u8, self.dom[c], UNK));
        self.dom[c] = 1 << v;
        self.val[c] = v;
        self.queue.push(c as u8);
        true
    }
    #[inline]
    fn unify(&mut self, c1: usize, c2: usize) -> bool {
        if c1 == c2 {
            return true;
        }
        let (v1, v2) = (self.val[c1], self.val[c2]);
        if v1 != UNK {
            return self.assign(c2, v1);
        }
        if v2 != UNK {
            return self.assign(c1, v2);
        }
        let d = self.dom[c1] & self.dom[c2];
        if d == 0 {
            return false;
        }
        if d & (d - 1) == 0 {
            let v = d.trailing_zeros() as u8;
            return self.assign(c1, v) && self.assign(c2, v);
        }
        if self.dom[c1] != d {
            self.set_dom(c1, d);
        }
        if self.dom[c2] != d {
            self.set_dom(c2, d);
        }
        true
    }
    fn undo(&mut self, mark: usize) {
        while self.trail.len() > mark {
            let (c, d, v) = self.trail.pop().unwrap();
            self.dom[c as usize] = d;
            self.val[c as usize] = v;
        }
    }
    fn propagate(&mut self, cx: &Ctx) -> bool {
        let n = cx.n;
        while let Some(c) = self.queue.pop() {
            let c = c as usize;
            if c >= n * n {
                continue; // depth cell
            }
            let (x, y) = (c / n, c % n);
            let v = self.val[c] as usize;
            if cx.anti {
                // T[pi y][pi x] = pi(T[x][y])
                let pc = cx.pi[y] as usize * n + cx.pi[x] as usize;
                if !self.assign(pc, cx.pi[v]) {
                    self.queue.clear();
                    return false;
                }
            }
            for t in 0..n {
                // (x y) t = x (y t)
                let q = self.val[y * n + t];
                if q != UNK && !self.unify(v * n + t, x * n + q as usize) {
                    self.queue.clear();
                    return false;
                }
                // (t x) y = t (x y)
                let p = self.val[t * n + x];
                if p != UNK && !self.unify(p as usize * n + y, t * n + v) {
                    self.queue.clear();
                    return false;
                }
            }
            for a in 0..n {
                for b in 0..n {
                    let w = self.val[a * n + b] as usize;
                    // a b = x:  (a b) y = a (b y)
                    if w == x {
                        let q = self.val[b * n + y];
                        if q != UNK && !self.assign(a * n + q as usize, v as u8) {
                            self.queue.clear();
                            return false;
                        }
                    }
                    // a b = y:  x (a b) = (x a) b
                    if w == y {
                        let p = self.val[x * n + a];
                        if p != UNK && !self.assign(p as usize * n + b, v as u8) {
                            self.queue.clear();
                            return false;
                        }
                    }
                }
            }
        }
        true
    }
}

type Hist = HashMap<(u8, u8, u32), u64>;

struct Run<'a> {
    cx: &'a Ctx,
    st: St,
    hist: Hist,
    task: u64,
    leaves: u64,
    nodes: u64,
    nosplit: bool,
    order: Vec<u8>,   // current cell order (the tail is re-sorted by depth sum in the nilpotent branch)
    csum: [u8; NC],   // depth(a) + depth(b) of cell (a, b) in the nilpotent branch
    nilmode: bool,
    cat: [u64; 3],
    exact: [u16; 16], // exact[k]: elements of depth k
    need: u16,        // elements of depth >= 2 (not the zero)
}

impl<'a> Run<'a> {
    /// Filter the surviving group elements against the current partial table.
    /// None: some image is smaller, so this node is not canonical.
    fn filter(&self, surv: &[(u32, u8)]) -> Option<Vec<(u32, u8)>> {
        let cx = self.cx;
        let n = cx.n;
        let nc = self.order.len();
        let val = &self.st.val;
        let mut out = Vec::with_capacity(surv.len());
        'g: for &(gi, ptr0) in surv {
            let g = &cx.g[gi as usize];
            let gi_ = &cx.ginv[gi as usize];
            let mut ptr = ptr0 as usize;
            while ptr < nc {
                let c = self.order[ptr] as usize;
                let tv = val[c];
                if tv == UNK {
                    break;
                }
                let (src, depth) = if c >= n * n {
                    (n * n + gi_[c - n * n] as usize, true)
                } else {
                    (gi_[c / n] as usize * n + gi_[c % n] as usize, false)
                };
                let sv = val[src];
                if sv == UNK {
                    break;
                }
                let im = if depth { sv } else { g[sv as usize] };
                if im < tv {
                    return None;
                }
                if im > tv {
                    continue 'g;
                }
                ptr += 1;
            }
            out.push((gi, ptr as u8));
        }
        Some(out)
    }

    /// Can the current partial table still be completed to something that is not 3-nilpotent
    /// (S^3 = {0})?  Only used when the diagonal is x^2 = z for... see `no3nil`.
    fn leaf(&mut self, stab: usize) {
        let n = self.cx.n;
        let v = &self.st.val;
        if self.cx.no3nil && (0..n).any(|t| v[n * n + t] == INF) {
            // nilpotent case: every element of depth >= 2 needs a factorization of exactly that depth
            for p in 0..n {
                let dp = v[n * n + p];
                if dp >= 2 && dp != INF {
                    let mut ok = false;
                    'f: for a in 0..n {
                        for b in 0..n {
                            if v[a * n + b] as usize == p && v[n * n + a] + v[n * n + b] == dp {
                                ok = true;
                                break 'f;
                            }
                        }
                    }
                    if !ok {
                        return;
                    }
                }
            }
        }
        let k = (0..n).filter(|&x| v[x * n + x] as usize == x).count() as u8;
        // 3-nilpotent or less: a zero z with (ab)c = z for all a, b, c
        let mut nil = 0u8;
        if let Some(z) = (0..n).find(|&z| (0..n).all(|t| v[z * n + t] as usize == z && v[t * n + z] as usize == z)) {
            nil = 1;
            'o: for a in 0..n {
                for b in 0..n {
                    let p = v[a * n + b] as usize;
                    for c in 0..n {
                        if v[p * n + c] as usize != z {
                            nil = 0;
                            break 'o;
                        }
                    }
                }
            }
        }
        // flags: 1 = S^3 is {0}, 2 = commutative, 4 = inverse (regular with commuting idempotents)
        let com = (0..n).all(|a| (0..a).all(|b| v[a * n + b] == v[b * n + a]));
        let regular = (0..n).all(|x| (0..n).any(|y| v[v[x * n + y] as usize * n + x] as usize == x));
        let idem: Vec<usize> = (0..n).filter(|&x| v[x * n + x] as usize == x).collect();
        let inv = regular && idem.iter().all(|&e| idem.iter().all(|&f| v[e * n + f] == v[f * n + e]));
        let flags = nil | (com as u8) << 1 | (inv as u8) << 2;
        *self.hist.entry((k, flags, stab as u32)).or_insert(0) += 1;
        self.leaves += 1;
        if !self.cx.anti && std::env::var("SD_LIST").is_ok() {
            // iso mode: print the canonical table (row-major digits)
            let t: String = (0..n * n).map(|c| (b'0' + v[c]) as char).collect();
            println!("table {t}");
        }
    }

    fn search(&mut self, depth: usize, start: usize, surv: &[(u32, u8)]) {
        let cx = self.cx;
        self.nodes += 1;
        if self.nilmode {
            self.cat[0] += 1;
        } else if (0..cx.n).any(|z| self.st.val[cx.n * cx.n + z] == KIND_GROUP) {
            self.cat[1] += 1;
        } else if self.st.val[cx.n * cx.n + cx.n - 1] == UNK && self.st.val[(cx.n - 1) * (cx.n + 1)] != UNK {
            self.cat[2] += 1;
        }
        if self.nilmode {
            // S^k = S . S^(k-1) = S^(k-1) . S: every element p of depth k >= 2 is (generator)(depth k-1) and
            // (depth k-1)(generator).  A filled cell gives that; otherwise probe the open cells that could.
            let n = cx.n;
            for side in 0..2 {
                let mut have: u16 = 0;
                for c in 0..n * n {
                    let v = self.st.val[c];
                    if v != UNK
                        && self.st.val[n * n + if side == 0 { c / n } else { c % n }] == 1
                        && self.csum[c] == self.st.val[n * n + v as usize]
                    {
                        have |= 1 << v;
                    }
                }
                let mut miss = self.need & !have;
                while miss != 0 {
                    let p = miss.trailing_zeros() as u8;
                    miss &= miss - 1;
                    let dpv = self.st.val[n * n + p as usize];
                    let mut found = false;
                    for c in 0..n * n {
                        if self.st.val[c] != UNK {
                            if self.st.val[c] == p
                                && self.csum[c] == dpv
                                && self.st.val[n * n + if side == 0 { c / n } else { c % n }] == 1
                            {
                                found = true; // filled by an earlier probe's consequence
                                break;
                            }
                            continue;
                        }
                        if self.csum[c] != dpv
                            || self.st.dom[c] & (1 << p) == 0
                            || self.st.val[n * n + if side == 0 { c / n } else { c % n }] != 1
                        {
                            continue;
                        }
                        let mark = self.st.trail.len();
                        let ok = self.st.assign(c, p) && self.st.propagate(cx);
                        self.st.queue.clear();
                        self.st.undo(mark);
                        if ok {
                            found = true;
                            break;
                        }
                        // p is impossible here
                        let nd = self.st.dom[c] & !(1 << p);
                        if nd == 0 {
                            return;
                        }
                        if nd & (nd - 1) == 0 {
                            if !(self.st.assign(c, nd.trailing_zeros() as u8) && self.st.propagate(cx)) {
                                self.st.queue.clear();
                                return;
                            }
                        } else {
                            self.st.set_dom(c, nd);
                        }
                    }
                    if !found {
                        return;
                    }
                }
            }
        }
        let nc = self.order.len();
        let mut pos = start;
        while pos < nc && self.st.val[self.order[pos] as usize] != UNK {
            pos += 1;
        }
        if cx.no3nil {
            // unique idempotent z marked "not a zero": its row must not be constant z
            let n = cx.n;
            if let Some(z) = (0..n).find(|&z| self.st.val[n * n + z] == KIND_GROUP) {
                if (0..n).all(|t| self.st.val[z * n + t] as usize == z) {
                    return;
                }
            }
        }
        if pos == nc {
            // a leaf above the split level is seen by every thread: thread 0 counts it
            if depth > cx.split || cx.tid == 0 {
                self.leaf(surv.len());
            }
            return;
        }
        let nosplit = std::mem::replace(&mut self.nosplit, false);
        if depth == cx.split && !nosplit {
            let t = self.task;
            self.task += 1;
            if t % cx.threads as u64 != cx.tid as u64 {
                return;
            }
        }
        let c = self.order[pos] as usize;
        if c >= cx.n * cx.n {
            self.depth_cells(depth, pos, c, surv);
            return;
        }
        let mut d = self.st.dom[c];
        while d != 0 {
            let v = d.trailing_zeros() as u8;
            d &= d - 1;
            let mark = self.st.trail.len();
            if self.st.assign(c, v) && self.st.propagate(cx) {
                if surv.len() == 1 {
                    self.search(depth + 1, pos + 1, surv); // only the identity is left
                } else if let Some(ns) = self.filter(surv) {
                    self.search(depth + 1, pos + 1, &ns);
                }
            }
            self.st.undo(mark);
        }
    }

    /// Depth cells.  If the diagonal says x^4 = z for the unique idempotent z, the semigroup is nilpotent
    /// and we branch on the depth of every element (longest factorization; INF for z).  Then
    /// a b = p != z forces depth(p) >= depth(a) + depth(b), and S^3 != {z} means some depth is >= 3.
    /// Otherwise all depth cells are 0.
    fn depth_cells(&mut self, depth: usize, pos: usize, c: usize, surv: &[(u32, u8)]) {
        let cx = self.cx;
        let n = cx.n;
        let x = c - n * n;
        let fixed: Vec<usize> = (0..n).filter(|&t| self.st.val[t * n + t] as usize == t).collect();
        let ntype = cx.no3nil
            && fixed.len() == 1
            && (0..n).all(|t| self.st.val[self.st.val[t * n + t] as usize * (n + 1)] as usize == fixed[0]);
        let mark = self.st.trail.len();
        if !ntype {
            for t in x..n {
                let ok = self.st.assign(n * n + t, 0);
                debug_assert!(ok);
            }
            self.st.queue.clear();
            self.nosplit = true;
            self.search(depth, pos, surv);
            self.st.undo(mark);
            return;
        }
        let z = fixed[0];
        if self.st.val[n * n + z] == UNK {
            // (B) z is not a zero: a nontrivial group, so not nilpotent.  Depth cells: KIND_GROUP at z, 0 elsewhere.
            for t in 0..n {
                let ok = self.st.assign(n * n + t, if t == z { KIND_GROUP } else { 0 });
                debug_assert!(ok);
            }
            self.st.queue.clear();
            self.nosplit = true;
            self.search(depth, pos, surv);
            self.st.undo(mark);
            // (A) z is a zero: nilpotent, branch on depths
            self.st.assign(n * n + z, INF);
            self.st.queue.clear();
            if x == z {
                self.after_depth(depth, pos, surv, true);
                self.st.undo(mark);
                return;
            }
        }
        for dv in 1..n as u8 {
            let mark = self.st.trail.len();
            let mut ok = self.st.assign(c, dv);
            // depth is invariant under the anti-automorphism
            if ok && cx.anti {
                let mut y = cx.pi[x] as usize;
                while ok && y != x {
                    ok = self.st.assign(n * n + y, dv);
                    y = cx.pi[y] as usize;
                }
            }
            self.st.queue.clear();
            if ok {
                // squares: depth(x^2) >= 2 depth(x) where both are known
                let v = &self.st.val;
                for t in 0..n {
                    let s = v[t * n + t] as usize;
                    let (dt, ds) = (v[n * n + t], v[n * n + s]);
                    if s != z && dt != UNK && ds != UNK && ds < 2 * dt {
                        ok = false;
                    }
                }
            }
            if ok {
                self.after_depth(depth + 1, pos, surv, false);
            }
            self.st.undo(mark);
        }
        self.st.undo(mark);
    }

    fn after_depth(&mut self, depth: usize, pos: usize, surv: &[(u32, u8)], same: bool) {
        let cx = self.cx;
        let n = cx.n;
        if (0..n).any(|t| self.st.val[n * n + t] == UNK) {
            if let Some(ns) = self.filter(surv) {
                self.nosplit = same;
                self.search(depth, pos, &ns);
            }
            return;
        }
        // all depths known: S^3 != {z}, and restrict every cell
        let dp: Vec<u8> = (0..n).map(|t| self.st.val[n * n + t]).collect();
        if !dp.iter().any(|&d| d >= 3 && d != INF) {
            return;
        }
        // a prefix of a longest factorization is a longest factorization: the depths are exactly 1..=m
        let m = dp.iter().filter(|&&d| d != INF).max().copied().unwrap();
        if (1..=m).any(|k| !dp.contains(&k)) {
            return;
        }
        let z = dp.iter().position(|&d| d == INF).unwrap();
        let mark = self.st.trail.len();
        let mut ok = true;
        'r: for a in 0..n {
            for b in 0..n {
                let need = if dp[a] == INF || dp[b] == INF { 200 } else { dp[a] + dp[b] };
                let mut mask: u16 = 1 << z;
                for p in 0..n {
                    if dp[p] != INF && dp[p] >= need {
                        mask |= 1 << p;
                    }
                }
                let c = a * n + b;
                if self.st.val[c] != UNK {
                    if mask & (1 << self.st.val[c]) == 0 {
                        ok = false;
                        break 'r;
                    }
                    continue;
                }
                let nd = self.st.dom[c] & mask;
                if nd == 0 {
                    ok = false;
                    break 'r;
                }
                if nd & (nd - 1) == 0 {
                    if !self.st.assign(c, nd.trailing_zeros() as u8) {
                        ok = false;
                        break 'r;
                    }
                } else if nd != self.st.dom[c] {
                    self.st.set_dom(c, nd);
                }
            }
        }
        if ok && self.st.propagate(cx) {
            // from here on compare and branch on the deepest products first
            let old = self.order.clone();
            let mut tail: Vec<u8> = old[2 * n..].to_vec();
            for &c in &tail {
                let (a, b) = (c as usize / n, c as usize % n);
                self.csum[c as usize] = if dp[a] == INF || dp[b] == INF { 255 } else { dp[a] + dp[b] };
            }
            for t in 0..n {
                self.csum[t * n + t] = if dp[t] == INF { 255 } else { 2 * dp[t] };
            }
            let cs = self.csum;
            tail.sort_by_key(|&c| (255 - cs[c as usize], c));
            self.order.truncate(2 * n);
            self.order.extend(tail);
            self.nilmode = true;
            self.exact = [0; 16];
            self.need = 0;
            for t in 0..n {
                if dp[t] != INF {
                    self.exact[dp[t] as usize] |= 1 << t;
                    if dp[t] >= 2 {
                        self.need |= 1 << t;
                    }
                }
            }
            let reset: Vec<(u32, u8)> = surv.iter().map(|&(g, p)| (g, p.min(2 * n as u8))).collect();
            if let Some(ns) = self.filter(&reset) {
                self.nosplit = same;
                let (n0, l0) = (self.nodes, self.leaves);
                self.search(depth, 2 * n, &ns);
                if std::env::var("SD_DEBUG").is_ok() && self.nodes - n0 > 0 {
                    let dg: Vec<u8> = (0..n).map(|t| self.st.val[t * n + t]).collect();
                    eprintln!("  depth {:?} diag {:?} surv {} nodes {} leaves {}", dp, dg, ns.len(), self.nodes - n0, self.leaves - l0);
                }
            }
            self.nilmode = false;
            self.order = old;
        } else {
            self.st.queue.clear();
        }
        self.st.undo(mark);
    }
}

fn perms(n: usize) -> Vec<[u8; MAXN]> {
    let mut out = vec![];
    let mut p: Vec<u8> = (0..n as u8).collect();
    fn rec(k: usize, p: &mut Vec<u8>, out: &mut Vec<[u8; MAXN]>) {
        let n = p.len();
        if k == n {
            let mut a = [0u8; MAXN];
            a[..n].copy_from_slice(p);
            out.push(a);
            return;
        }
        for i in k..n {
            p.swap(k, i);
            rec(k + 1, p, out);
            p.swap(k, i);
        }
    }
    rec(0, &mut p, &mut out);
    out
}

fn partitions(n: usize, max: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if n == 0 {
        out.push(cur.clone());
        return;
    }
    for k in (1..=n.min(max)).rev() {
        cur.push(k);
        partitions(n - k, k, cur, out);
        cur.pop();
    }
}

fn perm_of_type(lam: &[usize]) -> [u8; MAXN] {
    // fixed points get the smallest labels: cycles in increasing length
    let mut p = [0u8; MAXN];
    let mut s = 0;
    let mut l: Vec<usize> = lam.to_vec();
    l.sort();
    for &c in &l {
        for i in 0..c {
            p[s + i] = (s + (i + 1) % c) as u8;
        }
        s += c;
    }
    p
}

fn run_class(n: usize, anti: bool, pi: [u8; MAXN], threads: usize, split: usize, no3nil: bool) -> (Hist, u64, u64, usize) {
    let all = perms(n);
    let g: Vec<[u8; MAXN]> = if anti {
        all.into_iter().filter(|q| (0..n).all(|x| q[pi[x] as usize] == pi[q[x] as usize])).collect()
    } else {
        all
    };
    let ginv: Vec<[u8; MAXN]> = g
        .iter()
        .map(|q| {
            let mut r = [0u8; MAXN];
            for x in 0..n {
                r[q[x] as usize] = x as u8;
            }
            r
        })
        .collect();
    let mut order: Vec<u8> = (0..n).map(|x| (x * n + x) as u8).collect();
    order.extend((0..n).map(|x| (n * n + x) as u8));
    for x in 0..n {
        for y in 0..n {
            if x != y {
                order.push((x * n + y) as u8);
            }
        }
    }
    let gsize = g.len();
    let base = Ctx { n, anti, pi, order, g, ginv, threads, tid: 0, split, no3nil };
    let results: Vec<(Hist, u64, u64)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads)
            .map(|tid| {
                let mut cx = base.clone();
                cx.tid = tid;
                sc.spawn(move || {
                    let full: u16 = (1u16 << n) - 1;
                    let mut st = St { val: [UNK; NC], dom: [0; NC], trail: vec![], queue: vec![] };
                    for c in 0..n * n {
                        st.dom[c] = full;
                    }
                    for c in n * n..n * n + n {
                        st.dom[c] = 0xffff;
                    }
                    let surv: Vec<(u32, u8)> = (0..cx.g.len() as u32).map(|i| (i, 0u8)).collect();
                    let mut r = Run { cx: &cx, st, hist: HashMap::new(), task: 0, leaves: 0, nodes: 0, nosplit: false, order: cx.order.clone(), csum: [0; NC], nilmode: false, cat: [0; 3], exact: [0; 16], need: 0 };
                    r.search(0, 0, &surv);
                    eprintln!("  tid {} nodes: nil {} group {} depthvec {}", cx.tid, r.cat[0], r.cat[1], r.cat[2]);
                    (r.hist, r.leaves, r.nodes)
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut hist: Hist = HashMap::new();
    let (mut leaves, mut nodes) = (0, 0);
    for (h, l, nd) in results {
        for (k, v) in h {
            *hist.entry(k).or_insert(0) += v;
        }
        leaves += l;
        nodes += nd;
    }
    (hist, leaves, nodes, gsize)
}

#[cfg(windows)]
fn low_priority() {
    #[repr(C)]
    struct Pts {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
        fn SetPriorityClass(h: isize, class: u32) -> i32;
    }
    unsafe {
        let p = Pts { version: 1, control_mask: 1, state_mask: 0 }; // opt out of EcoQoS throttling
        SetProcessInformation(GetCurrentProcess(), 4, &p as *const _ as *const _, 12);
        SetPriorityClass(GetCurrentProcess(), 0x4000); // BELOW_NORMAL
    }
}
#[cfg(not(windows))]
fn low_priority() {}

fn main() {
    low_priority();
    let a: Vec<String> = std::env::args().collect();
    let mode = a.get(1).map(|s| s.as_str()).unwrap_or("");
    let n: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let threads: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(6);
    assert!(n >= 1 && n <= MAXN);
    let split = if n >= 7 { n + 2 } else { 0 };
    let threads = if split == 0 { 1 } else { threads };
    let t0 = Instant::now();
    match mode {
        "iso" => {
            let (hist, leaves, nodes, gs) = run_class(n, false, [0; MAXN], threads, split, false);
            let mut tot = 0f64;
            for (&(_, _, s), &c) in &hist {
                assert!(s as usize == 1 || gs % s as usize == 0);
                tot += c as f64 / 1.0 * 0.0 + 0.0;
                let _ = s;
            }
            let _ = tot;
            // without an anti-automorphism each leaf is one isomorphism class
            let mut byk = vec![0u64; n + 1];
            for (&(k, _, _), &c) in &hist {
                byk[k as usize] += c;
            }
            println!("iso n={n}: {} classes; by idempotents {:?}  (leaves {leaves}, nodes {nodes}, {:.1}s)",
                     byk.iter().sum::<u64>(), &byk[1..], t0.elapsed().as_secs_f64());
        }
        "sd" | "sdx" => {
            // sdx: skip the semigroups with S^3 = {0} (3-nilpotent, zero, trivial), known by formula
            let mut lams = vec![];
            partitions(n, n, &mut vec![], &mut lams);
            let first: usize = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(0);
            let last: usize = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(lams.len() - 1);
            for (i, lam) in lams.iter().enumerate() {
                if i < first || i > last {
                    continue;
                }
                let t1 = Instant::now();
                let pi = perm_of_type(lam);
                let (hist, leaves, nodes, gs) = run_class(n, true, pi, threads, split, mode == "sdx");
                let l: Vec<String> = lam.iter().map(|x| x.to_string()).collect();
                let mut keys: Vec<_> = hist.iter().collect();
                keys.sort();
                for (&(k, nil, s), &c) in keys {
                    println!("class {i} {} k={k} nil3={} com={} inv={} stab={s} count={c}", l.join(","), nil & 1, nil >> 1 & 1, nil >> 2 & 1);
                }
                eprintln!("class {i} [{}] |C|={gs} leaves={leaves} nodes={nodes} {:.1}s", l.join(","), t1.elapsed().as_secs_f64());
            }
            eprintln!("total {:.1}s", t0.elapsed().as_secs_f64());
        }
        _ => eprintln!("usage: selfdual sd|iso <n> [threads] [first_class] [last_class]"),
    }
}
