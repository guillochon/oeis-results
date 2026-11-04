//! A000112 (unlabeled posets) by Burnside over S_N:
//!
//!   N! a(N) = sum over cycle types lam of |C_lam| fix(lam).
//!
//! A g-invariant order makes every cycle of g an antichain, and the cycles carry a quotient order Q
//! (a labeled poset on the r cycles). For cycles u < v of lengths l_u, l_v, the invariant relation
//! between them is a nonempty set S_uv of residues mod d_uv = gcd(l_u, l_v): x_u(s) < x_v(t) iff
//! (t - s) mod d_uv lies in S_uv. Transitivity through u < w < v says: with g = gcd(d_uw, d_wv)
//! (which divides d_uv), every c mod d_uv with c mod g in (S_uw + S_wv) mod g must lie in S_uv.
//! So fix(g) = sum over labeled Q of the number W(Q) of such decorations. Only pairs with d > 1
//! have a choice, so if the cycle lengths are pairwise coprime, fix(g) = P(r) = A001035(r).
//!
//! Usage:
//!   a000112 burnside NMAX    a(1..NMAX) from the posets in ../../A006455/data (r <= 10)
//!                            and the larger-r values in ../stats.txt
//!   a000112 stat R L1,L2,..  fix for quotient size R and interacting cycle lengths L, from the
//!                            R-point posets on stdin (genposetg R t); prints a stats.txt line

use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

mod blocks;

const MAXN: usize = 16;

/// A001035(0..19): labeled posets.
const A001035: [u128; 20] = [
    1,
    1,
    3,
    19,
    219,
    4231,
    130023,
    6129859,
    431723379,
    44511042511,
    6611065248783,
    1396281677105899,
    414864951055853499,
    171850728381587059351,
    98484324257128207032183,
    77567171020440688353049939,
    83480529785490157813844256579,
    122152541250295322862941281269151,
    241939392597201176602897820148085023,
    0, // 646099441937791106493755218560442089979 does not fit in u128
];

/// A000112(0..16), for checks.
const A000112: [u128; 17] = [
    1,
    1,
    2,
    5,
    16,
    63,
    318,
    2045,
    16999,
    183231,
    2567284,
    46749427,
    1104891746,
    33823827452,
    1338193159771,
    68275077901156,
    4483130665195087,
];

/// A poset in natural labeling: down[j] is the strict down-set of j.
#[derive(Clone)]
struct Poset {
    n: usize,
    down: [u16; MAXN],
    cov: [u16; MAXN],
}

/// Parse one digraph6 line from genposetg's `t` output (Hasse diagram, edges x->y with x<y).
fn parse_d6(line: &str) -> Poset {
    let b = line.trim().as_bytes();
    assert!(b[0] == b'&', "not digraph6: {line}");
    let n = (b[1] - 63) as usize;
    assert!(n <= MAXN);
    let bit = |p: usize| -> bool {
        let c = (b[2 + p / 6] - 63) as u32;
        c >> (5 - p % 6) & 1 == 1
    };
    let mut cov = [0u16; MAXN];
    for x in 0..n {
        for y in 0..n {
            if bit(x * n + y) {
                assert!(x < y, "edge {x}->{y} not topological");
                cov[y] |= 1 << x;
            }
        }
    }
    let mut down = [0u16; MAXN];
    for y in 0..n {
        let mut d = 0u16;
        for x in 0..y {
            if cov[y] >> x & 1 == 1 {
                d |= down[x] | 1 << x;
            }
        }
        down[y] = d;
    }
    Poset { n, down, cov }
}

fn load(n: usize) -> Vec<Poset> {
    if n == 0 {
        return vec![Poset { n: 0, down: [0; MAXN], cov: [0; MAXN] }];
    }
    let path = format!("{}/../../A006455/data/posets{n}.d6", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines().filter(|l| !l.is_empty()).map(parse_d6).collect()
}

/// |Aut p| by backtracking over relation-preserving bijections.
fn automorphisms(p: &Poset) -> u128 {
    let n = p.n;
    let mut up = [0u16; MAXN];
    for y in 0..n {
        for x in 0..n {
            if p.down[y] >> x & 1 == 1 {
                up[x] |= 1 << y;
            }
        }
    }
    // Incomparable twins (same down- and up-sets) permute freely: |Aut| = Π k! · |Aut| of the
    // quotient, whose points carry their class size as a color.
    let mut rep: Vec<usize> = vec![];
    let mut size: Vec<u32> = vec![];
    let mut tw: u128 = 1;
    for x in 0..n {
        if let Some(c) = rep.iter().position(|&y| p.down[y] == p.down[x] && up[y] == up[x]) {
            size[c] += 1;
            tw *= size[c] as u128;
        } else {
            rep.push(x);
            size.push(1);
        }
    }
    let k = rep.len();
    let sig: Vec<(u32, u32, u32, u32)> = (0..k)
        .map(|i| { let x = rep[i]; (p.down[x].count_ones(), up[x].count_ones(), p.cov[x].count_ones(), size[i]) })
        .collect();
    fn rec(i: usize, k: usize, p: &Poset, rep: &[usize], sig: &[(u32, u32, u32, u32)], img: &mut [usize; MAXN], used: u16) -> u128 {
        if i == k {
            return 1;
        }
        let mut total = 0;
        for v in 0..k {
            if used >> v & 1 == 1 || sig[v] != sig[i] {
                continue;
            }
            let ok = (0..i).all(|j| {
                let (a, b, c, d) = (rep[i], rep[j], rep[v], rep[img[j]]);
                (p.down[a] >> b & 1) == (p.down[c] >> d & 1) && (p.down[b] >> a & 1) == (p.down[d] >> c & 1)
            });
            if ok {
                img[i] = v;
                total += rec(i + 1, k, p, rep, sig, img, used | 1 << v);
            }
        }
        total
    }
    let mut img = [0usize; MAXN];
    tw * rec(0, k, p, &rep, &sig, &mut img, 0)
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn fact(n: usize) -> u128 {
    (1..=n as u128).product()
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

fn class_size(lam: &[usize]) -> u128 {
    let n: usize = lam.iter().sum();
    let mut z: u128 = 1;
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for &c in lam {
        *counts.entry(c).or_default() += 1;
    }
    for (&c, &m) in &counts {
        z *= (c as u128).pow(m as u32) * fact(m);
    }
    fact(n) / z
}

/// The "interacting" cycles of lam: lengths > 1 sharing a factor with another such cycle.
fn interacting(lam: &[usize]) -> Vec<usize> {
    let nt: Vec<usize> = lam.iter().copied().filter(|&c| c > 1).collect();
    let mut t = vec![];
    for (i, &a) in nt.iter().enumerate() {
        if nt.iter().enumerate().any(|(j, &b)| j != i && gcd(a, b) > 1) {
            t.push(a);
        }
    }
    t
}

/// Mask of all residues mod d.
fn full(d: usize) -> u32 {
    (1u32 << d) - 1
}

/// Residues c mod d with c mod g in the set (A + B) mod g, where A is a set mod da, B mod db and
/// g divides da, db and d.
fn forced(a: u32, da: usize, b: u32, db: usize, g: usize, d: usize) -> u32 {
    let mut sum_mod_g = 0u32;
    for x in 0..da {
        if a >> x & 1 == 0 {
            continue;
        }
        for y in 0..db {
            if b >> y & 1 == 1 {
                sum_mod_g |= 1 << ((x + y) % g);
            }
        }
    }
    let mut res = 0u32;
    for c in 0..d {
        if sum_mod_g >> (c % g) & 1 == 1 {
            res |= 1 << c;
        }
    }
    res
}

/// Number of valid decorations W(Q, phi): marked points pts[i] carry cycle lengths len[i]; every other
/// point of q is coprime to all of them (it only forces).
fn weight(q: &Poset, pts: &[usize], len: &[usize]) -> u64 {
    let t = pts.len();
    let lt = |a: usize, b: usize| q.down[b] >> a & 1 == 1;
    let marked: u16 = pts.iter().fold(0, |m, &p| m | 1 << p);
    // Pair variables: (i, j) with pts[i] < pts[j] and d > 1.
    let mut var = vec![vec![usize::MAX; t]; t];
    let mut vars: Vec<(usize, usize, usize, bool)> = vec![]; // (i, j, d, forced_full)
    for i in 0..t {
        for j in 0..t {
            if i == j || !lt(pts[i], pts[j]) {
                continue;
            }
            let d = gcd(len[i], len[j]);
            if d == 1 {
                continue;
            }
            // Forced full by an unmarked point between, or a marked one with gcd(d_iw, d_wj) = 1.
            let between = q.down[pts[j]] & !q.down[pts[i]] & !(1 << pts[i]);
            let between = between & {
                // elements above pts[i]
                let mut up = 0u16;
                for x in 0..q.n {
                    if lt(pts[i], x) {
                        up |= 1 << x;
                    }
                }
                up
            };
            let mut ff = between & !marked != 0;
            for w in 0..t {
                if between >> pts[w] & 1 == 1 && gcd(gcd(len[i], len[w]), gcd(len[w], len[j])) == 1 {
                    ff = true;
                }
            }
            var[i][j] = vars.len();
            vars.push((i, j, d, ff));
        }
    }
    // Triples (i, w, j) of marked points with pts[i] < pts[w] < pts[j] where all three pairs carry
    // variables and g = gcd(d_iw, d_wj) > 1 (others either force full or are trivial).
    let mut triples: Vec<(usize, usize, usize, usize)> = vec![]; // (var iw, var wj, var ij, g)
    for &(i, j, _, ff) in &vars {
        if ff {
            continue;
        }
        for w in 0..t {
            if var[i][w] != usize::MAX && var[w][j] != usize::MAX {
                let g = gcd(vars[var[i][w]].2, vars[var[w][j]].2);
                if g > 1 {
                    triples.push((var[i][w], var[w][j], var[i][j], g));
                }
            }
        }
    }
    let nv = vars.len();
    // Check each triple as soon as its last variable is assigned.
    let mut at: Vec<Vec<(usize, usize, usize, usize)>> = vec![vec![]; nv];
    for &tr in &triples {
        at[tr.0.max(tr.1).max(tr.2)].push(tr);
    }
    let mut val = vec![0u32; nv];
    fn rec(k: usize, vars: &[(usize, usize, usize, bool)], at: &[Vec<(usize, usize, usize, usize)>], val: &mut Vec<u32>) -> u64 {
        if k == vars.len() {
            return 1;
        }
        let d = vars[k].2;
        let (lo, hi) = if vars[k].3 { (full(d), full(d)) } else { (1, full(d)) };
        let mut s = 0;
        'm: for m in lo..=hi {
            val[k] = m;
            for &(a, b, c, g) in &at[k] {
                if forced(val[a], vars[a].2, val[b], vars[b].2, g, vars[c].2) & !val[c] != 0 {
                    continue 'm;
                }
            }
            s += rec(k + 1, vars, at, val);
        }
        s
    }
    rec(0, &vars, &at, &mut val)
}

/// Key for the W cache: for each ordered pair of marked points, 0 = not below,
/// 1 = below with no forcing point between (unmarked, or marked with gcd 1), 2 = below and forced.
fn pattern(q: &Poset, pts: &[usize], up: &[u16; MAXN], cop: &[[u8; 8]; 8]) -> u64 {
    let t = pts.len();
    let marked: u16 = pts.iter().fold(0, |m, &p| m | 1 << p);
    let mut key = 0u64;
    for i in 0..t {
        for j in 0..t {
            let mut c = 0;
            if i != j && q.down[pts[j]] >> pts[i] & 1 == 1 {
                c = 1;
                let between = q.down[pts[j]] & up[pts[i]];
                if between & !marked != 0 {
                    c = 2;
                } else {
                    let mut w = cop[i][j];
                    while w != 0 {
                        let k = w.trailing_zeros() as usize;
                        w &= w - 1;
                        if between >> pts[k] & 1 == 1 {
                            c = 2;
                            break;
                        }
                    }
                }
            }
            key = key * 3 + c;
        }
    }
    key
}

/// Fast sum_phi for two or three 2-cycles (the only classes needed on the largest quotients).
/// Two cycles with d = gcd: W = 2^d - 1 on a cover, else 1, so the sum is r(r-1) + 2(2^d - 2)E
/// (E = Hasse edges).
/// Three: away from 3-chains W = 3^(covers in the triple); on a chain a < b < c whose interval (a, c)
/// is exactly {b}, W = 13 (S_ac must contain S_ab + S_bc) instead of 9. Hasse diagrams have no
/// triangles, so the sum over unordered triples is C(r,3) + 2E(r-2) + 4 sum_v C(deg v, 2) + 4 L,
/// with L the number of such lonely chains; ordered maps are 6 times that.
fn sum_phi_fast(q: &Poset, len: &[usize]) -> Option<u64> {
    let r = q.n as u64;
    let e: u64 = (0..q.n).map(|y| q.cov[y].count_ones() as u64).sum();
    if (len.len() == 4 || len.len() == 5) && len.iter().all(|&x| x == 2) {
        return Some(all2(q, len.len()));
    }
    match len {
        // two interacting cycles: W = 2^d - 1 on a cover, else 1
        &[a, b] => Some(r * (r - 1) + 2 * ((1u64 << gcd(a, b)) - 2) * e),
        [2, 2, 2] => {
            let mut deg = [0u64; MAXN];
            for y in 0..q.n {
                deg[y] += q.cov[y].count_ones() as u64;
                for x in 0..q.n {
                    if q.cov[y] >> x & 1 == 1 {
                        deg[x] += 1;
                    }
                }
            }
            let paths: u64 = deg[..q.n].iter().map(|&d| d * d.saturating_sub(1) / 2).sum();
            // lonely chains: a covered by b covered by c, and b the only element strictly between a, c
            let mut lonely = 0u64;
            for c in 0..q.n {
                let mut cb = q.cov[c];
                while cb != 0 {
                    let b = cb.trailing_zeros() as usize;
                    cb &= cb - 1;
                    let mut ba = q.cov[b];
                    while ba != 0 {
                        let a = ba.trailing_zeros() as usize;
                        ba &= ba - 1;
                        // elements x with a < x < c: below c and above a
                        let mut between = 0u32;
                        for x in 0..q.n {
                            if q.down[c] >> x & 1 == 1 && q.down[x] >> a & 1 == 1 {
                                between += 1;
                            }
                        }
                        if between == 1 {
                            lonely += 1;
                        }
                    }
                }
            }
            let tri = r * (r - 1) * (r - 2) / 6 + 2 * e * (r - 2) + 4 * paths + 4 * lonely;
            Some(6 * tri)
        }
        _ => None,
    }
}

thread_local! {
    static W_CACHE: std::cell::RefCell<HashMap<u64, u64>> = std::cell::RefCell::new(HashMap::new());
}

/// t two-cycles (t = 4, 5): t! Σ over t-sets M of W(M). All d = 2, so a cover pair has 3
/// decorations and W(M) = 3^e(M) unless M holds a "long" unforced pair (u < v, not a cover, the
/// interval (u,v) inside M). With exactly one long pair, a lonely chain a ⋖ b ⋖ c with (a,c) = {b},
/// W(M) = 13·3^(e(M)-2). So
///   Σ W = Σ_M 3^e(M) + Σ_{lonely L} Σ_{M ⊇ L} 4·3^(e(M)-2)
///         + Σ_{multi M} [W(M) − 3^e(M) − #lonely(M)·4·3^(e(M)-2)],
/// where multi M (two or more long pairs, or one with an interval of 2+ points) are found from
/// their cores and done exactly with `weight`.
fn all2(q: &Poset, t: usize) -> u64 {
    let n = q.n;
    let r = n as u64;
    let mut up = [0u16; MAXN];
    let mut adj = [0u16; MAXN];
    for y in 0..n {
        for x in 0..n {
            if q.down[y] >> x & 1 == 1 {
                up[x] |= 1 << y;
            }
            if q.cov[y] >> x & 1 == 1 {
                adj[x] |= 1 << y;
                adj[y] |= 1 << x;
            }
        }
    }
    let p3 = |k: u32| -> u64 { 3u64.pow(k) };
    let edges_in = |m: u16| -> u32 {
        let mut s = 0;
        let mut b = m;
        while b != 0 {
            let x = b.trailing_zeros() as usize;
            b &= b - 1;
            s += (adj[x] & m).count_ones();
        }
        s / 2
    };
    let c = |a: u64, b: u64| -> u64 { if b > a { 0 } else { (0..b).fold(1u64, |acc, i| acc * (a - i) / (i + 1)) } };
    // Σ_M 3^e(M)
    let term_a: u64 = if t == 4 {
        let deg: Vec<u64> = (0..n).map(|x| adj[x].count_ones() as u64).collect();
        let e: u64 = deg.iter().sum::<u64>() / 2;
        let p2: u64 = deg.iter().map(|&d| d * d.saturating_sub(1) / 2).sum();
        let s3: u64 = deg.iter().map(|&d| c(d, 3)).sum();
        let mut p4 = 0u64;
        for y in 0..n {
            let mut b = adj[y];
            while b != 0 {
                let x = b.trailing_zeros() as usize;
                b &= b - 1;
                if x < y {
                    p4 += (deg[x] - 1) * (deg[y] - 1);
                }
            }
        }
        let mut c4x2 = 0u64;
        for u in 0..n {
            for v in u + 1..n {
                c4x2 += c((adj[u] & adj[v]).count_ones() as u64, 2);
            }
        }
        c(r, 4) + 2 * e * c(r - 2, 2) + 4 * p2 * (r - 3) + 4 * (e * e.saturating_sub(1) / 2 - p2) + 8 * s3 + 8 * p4 + 16 * (c4x2 / 2)
    } else {
        let mut s = 0u64;
        let mut m: u32 = (1 << t) - 1;
        while m < 1 << n {
            s += p3(edges_in(m as u16));
            let cc = m & m.wrapping_neg();
            let rr = m + cc;
            m = (((rr ^ m) >> 2) / cc) | rr;
        }
        s
    };
    // lonely chains and long pairs with bigger intervals
    let mut lonely: Vec<u16> = vec![];
    let mut cores: Vec<u16> = vec![];
    for cc in 0..n {
        let mut b = q.down[cc] & !q.cov[cc];
        while b != 0 {
            let a = b.trailing_zeros() as usize;
            b &= b - 1;
            let iv = up[a] & q.down[cc];
            let k = iv.count_ones() as usize;
            let core = iv | 1 << a | 1 << cc;
            if k == 1 {
                lonely.push(core);
            } else if k + 2 <= t {
                cores.push(core);
            }
        }
    }
    let mut term_b = 0u64;
    let all: u16 = ((1u32 << n) - 1) as u16;
    for &l in &lonely {
        let rest = all & !l;
        if t == 4 {
            let mut b = rest;
            while b != 0 {
                let x = b.trailing_zeros() as usize;
                b &= b - 1;
                term_b += 4 * p3((adj[x] & l).count_ones());
            }
        } else {
            let xs: Vec<usize> = (0..n).filter(|&x| rest >> x & 1 == 1).collect();
            for i in 0..xs.len() {
                for j in i + 1..xs.len() {
                    let (x, y) = (xs[i], xs[j]);
                    let k = (adj[x] & l).count_ones() + (adj[y] & l).count_ones() + (adj[x] >> y & 1) as u32;
                    term_b += 4 * p3(k);
                }
            }
        }
    }
    for i in 0..lonely.len() {
        for j in i + 1..lonely.len() {
            let u = lonely[i] | lonely[j];
            if (u.count_ones() as usize) <= t {
                cores.push(u);
            }
        }
    }
    // multi sets: every t-set containing a core
    let mut multi: Vec<u16> = vec![];
    for &core in &cores {
        let need = t - core.count_ones() as usize;
        let rest: Vec<usize> = (0..n).filter(|&x| core >> x & 1 == 0).collect();
        let mut pick = vec![0usize; need];
        fn rec(k: usize, start: usize, rest: &[usize], pick: &mut Vec<usize>, core: u16, out: &mut Vec<u16>) {
            if k == pick.len() {
                out.push(pick.iter().fold(core, |m, &x| m | 1 << x));
                return;
            }
            for i in start..rest.len() {
                pick[k] = rest[i];
                rec(k + 1, i + 1, rest, pick, core, out);
            }
        }
        rec(0, 0, &rest, &mut pick, core, &mut multi);
    }
    multi.sort_unstable();
    multi.dedup();
    let mut corr: i64 = 0;
    let twos = vec![2usize; t];
    for &m in &multi {
        let pts: Vec<usize> = (0..n).filter(|&x| m >> x & 1 == 1).collect();
        // W depends only on the pattern of M (all outside points are fixed points)
        let key = pattern(q, &pts, &up, &[[0u8; 8]; 8]) * 8 + t as u64;
        let w = W_CACHE.with(|c| *c.borrow_mut().entry(key).or_insert_with(|| weight(q, &pts, &twos))) as i64;
        let e = edges_in(m);
        let nl = lonely.iter().filter(|&&l| l & m == l).count() as i64;
        corr += w - p3(e) as i64 - nl * 4 * if e >= 2 { p3(e - 2) as i64 } else { 0 };
    }
    let total = (term_a + term_b) as i64 + corr;
    let tf: u64 = (1..=t as u64).product();
    total as u64 * tf
}

/// Sum over injective maps of the marked cycles into q of W.
fn sum_phi(q: &Poset, len: &[usize], cache: &mut HashMap<u64, u64>) -> u64 {
    // up-sets and, per ordered marked pair (i, j), the marked w whose cycle forces (i, j) when between
    let mut up = [0u16; MAXN];
    for y in 0..q.n {
        for x in 0..q.n {
            if q.down[y] >> x & 1 == 1 {
                up[x] |= 1 << y;
            }
        }
    }
    let t = len.len();
    assert!(t <= 8);
    let mut cop = [[0u8; 8]; 8];
    for i in 0..t {
        for j in 0..t {
            for w in 0..t {
                if w != i && w != j && gcd(gcd(len[i], len[w]), gcd(len[w], len[j])) == 1 {
                    cop[i][j] |= 1 << w;
                }
            }
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn rec(i: usize, q: &Poset, len: &[usize], pts: &mut Vec<usize>, used: u16, cache: &mut HashMap<u64, u64>, up: &[u16; MAXN], cop: &[[u8; 8]; 8]) -> u64 {
        if i == len.len() {
            let key = pattern(q, pts, up, cop);
            return *cache.entry(key).or_insert_with(|| weight(q, pts, len));
        }
        // W is symmetric among marked cycles of equal length (len is sorted, so they are adjacent):
        // take those in increasing point order and multiply by the multiplicity factorials below.
        let start = if i > 0 && len[i] == len[i - 1] { pts[i - 1] + 1 } else { 0 };
        let mut s = 0;
        for p in start..q.n {
            if used >> p & 1 == 0 {
                pts.push(p);
                s += rec(i + 1, q, len, pts, used | 1 << p, cache, up, cop);
                pts.pop();
            }
        }
        s
    }
    let mut mult = 1u64;
    let mut run = 1u64;
    for i in 1..len.len() {
        run = if len[i] == len[i - 1] { run + 1 } else { 1 };
        mult *= run;
    }
    debug_assert!(len.windows(2).all(|w| w[0] >= w[1]));
    mult * rec(0, q, len, &mut vec![], 0, cache, &up, &cop)
}

fn burnside(nmax: usize) {
    disable_power_throttling();
    let t0 = Instant::now();
    let rmax = 10;
    let posets: Vec<Vec<(Poset, u128)>> = (0..=rmax)
        .map(|r| load(r).into_iter().map(|p| { let w = fact(r) / automorphisms(&p); (p, w) }).collect())
        .collect();
    for r in 0..=rmax {
        let s: u128 = posets[r].iter().map(|x| x.1).sum();
        assert_eq!(s, A001035[r], "labeled check r={r}");
    }
    eprintln!("loaded posets r<=10 with |Aut| ({:.1}s)", t0.elapsed().as_secs_f64());
    // Statistics for r > 10 computed by `stat` (streamed posets), one "r l1,l2,.. value" per line.
    let mut extra: BTreeMap<(usize, Vec<usize>), u128> = BTreeMap::new();
    if let Ok(text) = std::fs::read_to_string(format!("{}/../stats.txt", env!("CARGO_MANIFEST_DIR"))) {
        for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
            let f: Vec<&str> = line.split_whitespace().collect();
            let lens = f[1].split(',').map(|x| x.parse().unwrap()).collect();
            let key = (f[0].parse().unwrap(), lens);
            let v: u128 = f[2].parse().unwrap();
            if let Some(&old) = extra.get(&key) {
                assert_eq!(old, v, "stats.txt disagrees with itself on {key:?}");
            }
            extra.insert(key, v);
        }
    }
    // fix cache keyed by (r, interacting lengths), persisted in fix_small.txt ("r lens value" lines) so a
    // rerun skips the r <= 10 sums. A000112_TRUST_FAST=1 uses the closed forms without recomputing the
    // generic count (they were asserted equal on every r <= 10 class in the full run).
    let cache_path = format!("{}/../fix_small.txt", env!("CARGO_MANIFEST_DIR"));
    let mut cache: BTreeMap<(usize, Vec<usize>), u128> = BTreeMap::new();
    if let Ok(text) = std::fs::read_to_string(&cache_path) {
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let f: Vec<&str> = line.split_whitespace().collect();
            let lens = f[1].split(',').map(|x| x.parse().unwrap()).collect();
            cache.insert((f[0].parse().unwrap(), lens), f[2].parse().unwrap());
        }
    }
    let trust_fast = std::env::var("A000112_TRUST_FAST").is_ok();
    let threads = 6;
    for n in 1..=nmax {
        let mut lams = vec![];
        partitions(n, n, &mut vec![], &mut lams);
        let mut total: u128 = 0;
        let mut hard: u128 = 0;
        for lam in &lams {
            let r = lam.len();
            let t = interacting(lam);
            let fx = if t.is_empty() {
                A001035[r]
            } else if let Some(&v) = extra.get(&(r, t.clone())) {
                v
            } else {
                assert!(r <= rmax, "need posets on {r} points for {lam:?}");
                if let Some(&v) = cache.get(&(r, t.clone())) {
                    v
                } else {
                    // split the posets over threads, each with its own pattern cache
                    let ps = &posets[r];
                    let chunk = ps.len().div_ceil(threads).max(1);
                    let tt = &t;
                    let s: u128 = std::thread::scope(|sc| {
                        let hs: Vec<_> = ps
                            .chunks(chunk)
                            .map(|part| {
                                sc.spawn(move || {
                                    let mut wc = HashMap::new();
                                    part.iter()
                                        .map(|(q, w)| {
                                            let fast = sum_phi_fast(q, tt);
                                            let g = match fast {
                                                Some(f) if trust_fast => f,
                                                _ => {
                                                    let g = sum_phi(q, tt, &mut wc);
                                                    if let Some(f) = fast {
                                                        assert_eq!(f, g, "fast path {tt:?} r={r}");
                                                    }
                                                    g
                                                }
                                            };
                                            *w * g as u128
                                        })
                                        .sum::<u128>()
                                })
                            })
                            .collect();
                        hs.into_iter().map(|h| h.join().unwrap()).sum()
                    });
                    let falling: u128 = ((r - t.len() + 1)..=r).map(|x| x as u128).product();
                    assert!(s % falling == 0);
                    let v = s / falling;
                    cache.insert((r, t.clone()), v);
                    let l: Vec<String> = t.iter().map(|x| x.to_string()).collect();
                    use std::io::Write;
                    let mut fh = std::fs::OpenOptions::new().create(true).append(true).open(&cache_path).unwrap();
                    writeln!(fh, "{r} {} {v}", l.join(",")).unwrap();
                    v
                }
            };
            let c = class_size(lam) * fx;
            total += c;
            if !t.is_empty() {
                hard += c;
                // per-class record of the hard classes (cycle type, its fix), for the shares chart
                let l: Vec<String> = lam.iter().map(|x| x.to_string()).collect();
                println!("class {n} {} {fx}", l.join(","));
            }
        }
        let f = fact(n);
        assert!(total % f == 0, "n={n}: not divisible");
        let a = total / f;
        let ok = if n < A000112.len() { if a == A000112[n] { "OK" } else { "MISMATCH" } } else { "new?" };
        println!(
            "a({n}) = {a}  [{ok}]  hard classes {:.4}%  ({:.1}s)",
            100.0 * hard as f64 / total as f64,
            t0.elapsed().as_secs_f64()
        );
    }
}

/// fix for quotient size r and interacting lengths `lens`, from posets on r points read from stdin.
fn stat(r: usize, keys: Vec<Vec<usize>>) {
    use std::io::BufRead;
    disable_power_throttling();
    let t0 = Instant::now();
    let fr = fact(r);
    let mut wc: Vec<HashMap<u64, u64>> = keys.iter().map(|_| HashMap::new()).collect();
    let mut sums = vec![0u128; keys.len()];
    let (mut cnt, mut lab) = (0u64, 0u128);
    // Checkpoints (A000112_CKPT=file): genposetg's order is deterministic, so a restarted slice
    // skips the posets already summed and resumes from the saved totals.
    let ckpt = std::env::var("A000112_CKPT").ok();
    let mut skip = 0u64;
    if let Some(path) = &ckpt {
        if let Ok(text) = std::fs::read_to_string(path) {
            let f: Vec<&str> = text.split_whitespace().collect();
            assert_eq!(f.len(), 2 + keys.len(), "checkpoint {path} does not match the keys");
            skip = f[0].parse().unwrap();
            lab = f[1].parse().unwrap();
            for k in 0..keys.len() {
                sums[k] = f[2 + k].parse().unwrap();
            }
            eprintln!("resuming after {skip} posets");
        }
    }
    let every: u64 = std::env::var("A000112_CKPT_EVERY").ok().and_then(|v| v.parse().ok()).unwrap_or(50_000_000);
    let mut seen = 0u64;
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        if !line.starts_with('&') {
            continue;
        }
        seen += 1;
        if seen <= skip {
            continue;
        }
        let q = parse_d6(&line);
        assert_eq!(q.n, r);
        let w = fr / automorphisms(&q);
        lab += w;
        for (k, lens) in keys.iter().enumerate() {
            let v = match sum_phi_fast(&q, lens) { Some(f) => f, None => sum_phi(&q, lens, &mut wc[k]) };
            sums[k] += w * v as u128;
        }
        cnt += 1;
        if let Some(path) = &ckpt {
            if seen % every == 0 {
                let mut t = format!("{seen} {lab}");
                for v in &sums {
                    t += &format!(" {v}");
                }
                let tmp = format!("{path}.tmp");
                std::fs::write(&tmp, t).unwrap();
                std::fs::rename(&tmp, path).unwrap();
            }
        }
    }
    cnt += skip;
    eprintln!("stat {r}: {cnt} posets, {:.1}s", t0.elapsed().as_secs_f64());
    // raw sums, combined over slices by combine_stats.py: fix = sum / (r)_t
    println!("lab {r} {lab} {cnt}");
    for (k, lens) in keys.iter().enumerate() {
        let l: Vec<String> = lens.iter().map(|x| x.to_string()).collect();
        println!("raw {r} {} {}", l.join(","), sums[k]);
    }
}

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
        fn SetPriorityClass(h: isize, class: u32) -> i32;
    }
    let st = St { version: 1, control_mask: 1, state_mask: 0 };
    unsafe {
        SetProcessInformation(GetCurrentProcess(), 4, &st as *const _ as *const std::ffi::c_void, std::mem::size_of::<St>() as u32);
        // long sweeps run below normal priority (WSL-launched processes don't inherit it)
        SetPriorityClass(GetCurrentProcess(), 0x4000);
    }
}
#[cfg(not(windows))]
fn disable_power_throttling() {}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("agg") => blocks::cmd_agg(args[2].parse().unwrap()),
        Some("blocks") => blocks::cmd_blocks(args[2].parse().unwrap(), args.get(3).map_or(10, |s| s.parse().unwrap())),
        Some("blocks1") => blocks::cmd_blocks1(args[2].parse().unwrap(), &args[3], args.get(4).map_or(10, |s| s.parse().unwrap())),
        Some("stat") => stat(
            args[2].parse().unwrap(),
            args[3].split('/').map(|k| k.split(',').map(|x| x.parse().unwrap()).collect()).collect(),
        ),
        Some("burnside") => burnside(args.get(2).map_or(12, |s| s.parse().unwrap())),
        _ => eprintln!("usage: a000112 burnside NMAX"),
    }
}
