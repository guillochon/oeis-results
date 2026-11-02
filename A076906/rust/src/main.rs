//! Transfer-matrix count of meanders crossing two perpendicular roads
//! (A076906, A076907, A077551 and relatives).  Same algorithm as tm.py.
//!
//! The processed region lies left of a cut curve.  The cut is a sequence of
//! symbols, top to bottom: markers (where the cut crosses a ray of the axes) and
//! items (river edges crossing the cut).  Paired items connect inside the processed
//! region and are non-crossing, so they are written as parentheses; START / END
//! items lead to the river's ends at infinity.  Markers are anonymous; the stage
//! says which rays they are:
//!   stage 1 (A): [W]       W points, cut x = c < 0 moving right
//!   stage 2 (B): [N, W]    N points, N marker moving down from +inf
//!   stage 3 (C): [E, S]    S points, S marker moving down from O
//!   stage 4 (D): [E]       E points, cut x = c > 0 moving right
//! A point at the active marker has one half-edge on each side; each closes the
//! adjacent item or opens a new one.  Closed loops are forbidden.
//!
//! Usage: meander2 NMAX END [start=q2|q3] [first] [noW|noN|noS|noE] [endfirst]
//!        END in q1 q2 q3 q4
//! Prints, for each n <= NMAX, the count of paths starting at infinity in Q3 and
//! ending in END, and with `first`: the part whose first crossing is on the x axis,
//! and the part whose first crossings from both ends are on the x axis (END in Q2/Q3).

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!()
    }
    fn write_u128(&mut self, k: u128) {
        let h = (k as u64) ^ ((k >> 64) as u64).rotate_left(29);
        self.0 = h.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (h >> 31);
    }
}
type Map = HashMap<u128, u128, BuildHasherDefault<Fx>>;

const OPEN: u8 = 1;
const CLOSE: u8 = 2;
const START: u8 = 3;
const END: u8 = 4;
const MARK: u8 = 5;
const MAXLEN: usize = 40;
/// Track which axis the river crosses first (only for checking A077551, which by
/// reflection in y = x equals A076907(n)/2 for even n and A076906(n) for odd n).
static TRACK_FIRST: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
/// Bit `stage` set: the ray processed in that stage carries no points (used for the
/// reflection-symmetric half curves; see README).
static FORBID: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);

#[derive(Clone)]
struct State {
    stage: u8,
    done: bool,
    // axis of the first crossing seen from the START end (bits 0-1) and from the END
    // end (bits 2-3): 0 not yet, 1 x axis, 2 y axis
    first: u8,
    seq: Vec<u8>,
}

fn encode(s: &State) -> u128 {
    assert!(s.seq.len() <= MAXLEN, "cut too long");
    let mut k: u128 = s.stage as u128 | (s.done as u128) << 3 | (s.first as u128) << 4;
    const _: () = assert!(8 + 3 * MAXLEN <= 128);
    for (i, &x) in s.seq.iter().enumerate() {
        k |= (x as u128) << (8 + 3 * i);
    }
    k
}

fn decode(k: u128) -> State {
    let mut seq = Vec::with_capacity(MAXLEN);
    let mut i = 0;
    loop {
        let x = ((k >> (8 + 3 * i)) & 7) as u8;
        if x == 0 || i == MAXLEN {
            break;
        }
        seq.push(x);
        i += 1;
    }
    State { stage: (k & 7) as u8, done: (k >> 3) & 1 == 1, first: ((k >> 4) & 15) as u8, seq }
}

fn is_item(x: u8) -> bool {
    x != MARK
}

fn marks(seq: &[u8]) -> Vec<usize> {
    seq.iter().enumerate().filter(|&(_, &x)| x == MARK).map(|(i, _)| i).collect()
}

fn active(s: &State) -> usize {
    let m = marks(&s.seq);
    match s.stage {
        1 | 2 | 4 => m[0],
        3 => m[1],
        _ => unreachable!(),
    }
}

/// Partner of the parenthesis at j (markers and START/END are skipped).
fn partner(seq: &[u8], j: usize) -> usize {
    let mut d = 0i32;
    if seq[j] == OPEN {
        for k in j..seq.len() {
            match seq[k] {
                OPEN => d += 1,
                CLOSE => {
                    d -= 1;
                    if d == 0 {
                        return k;
                    }
                }
                _ => {}
            }
        }
    } else {
        for k in (0..=j).rev() {
            match seq[k] {
                CLOSE => d += 1,
                OPEN => {
                    d -= 1;
                    if d == 0 {
                        return k;
                    }
                }
                _ => {}
            }
        }
    }
    unreachable!()
}

/// Side x axis (1) or y axis (2) of the ray the active marker sits on.
fn ray_axis(stage: u8) -> u8 {
    match stage {
        1 | 4 => 1,
        _ => 2,
    }
}

fn step(s: &State, out: &mut Vec<State>) {
    if s.done || FORBID.load(std::sync::atomic::Ordering::Relaxed) >> s.stage & 1 == 1 {
        return;
    }
    let i = active(s);
    let seq = &s.seq;
    let before = i > 0 && is_item(seq[i - 1]);
    let after = i + 1 < seq.len() && is_item(seq[i + 1]);
    let mk = |seq: Vec<u8>, done: bool, first: u8| State { stage: s.stage, done, first, seq };
    let track = TRACK_FIRST.load(std::sync::atomic::Ordering::Relaxed);
    // record the axis when the START / END edge is closed for the first time
    let upd = |f: u8, x: u8| match x {
        START if track && f & 3 == 0 => f | ray_axis(s.stage),
        END if track && f >> 2 == 0 => f | ray_axis(s.stage) << 2,
        _ => f,
    };

    // open both
    let mut v = seq[..i].to_vec();
    v.extend_from_slice(&[OPEN, MARK, CLOSE]);
    v.extend_from_slice(&seq[i + 1..]);
    out.push(mk(v, false, s.first));
    // close one side, open the other: the edge passes through the point
    if before {
        let mut v = seq.clone();
        v.swap(i - 1, i);
        out.push(mk(v, false, upd(s.first, seq[i - 1])));
    }
    if after {
        let mut v = seq.clone();
        v.swap(i, i + 1);
        out.push(mk(v, false, upd(s.first, seq[i + 1])));
    }
    if before && after {
        let (a, b) = (seq[i - 1], seq[i + 1]);
        let first = upd(upd(s.first, a), b);
        let mut v = seq.clone();
        match (a, b) {
            (OPEN, CLOSE) => return, // they are partners: closed loop
            (CLOSE, OPEN) => {}      // outer partners stay ( ... )
            (OPEN, OPEN) => {
                let q = partner(seq, i + 1);
                v[q] = OPEN;
            }
            (CLOSE, CLOSE) => {
                let p = partner(seq, i - 1);
                v[p] = CLOSE;
            }
            (START, END) | (END, START) => {
                v.remove(i + 1);
                v.remove(i - 1);
                if v.iter().any(|&x| is_item(x)) {
                    return;
                }
                out.push(mk(v, true, first));
                return;
            }
            (lab, p) | (p, lab) if (lab == START || lab == END) && (p == OPEN || p == CLOSE) => {
                let j = if a == lab { i + 1 } else { i - 1 };
                let q = partner(seq, j);
                v[q] = lab;
            }
            _ => unreachable!(),
        }
        v.remove(i + 1);
        v.remove(i - 1);
        out.push(mk(v, false, first));
    }
}

/// Free stage changes (no point processed).
fn transition(s: &State) -> Option<State> {
    let m = marks(&s.seq);
    match s.stage {
        1 => {
            let mut v = vec![MARK];
            v.extend_from_slice(&s.seq);
            Some(State { stage: 2, seq: v, ..s.clone() })
        }
        2 if m[1] == m[0] + 1 => Some(State { stage: 3, ..s.clone() }),
        3 if m[1] == s.seq.len() - 1 => {
            let mut v = s.seq.clone();
            v.pop();
            Some(State { stage: 4, seq: v, ..s.clone() })
        }
        _ => None,
    }
}

/// Can the items still be closed with `rem` more points?  Each point closes at
/// most one item on each of its two sides; `spare` items may survive to infinity.
fn feasible(s: &State, rem: usize, spare: usize, endq: &str) -> bool {
    // Q1 items close at N or E points, Q2 at W or N, Q3 at W or S, Q4 at S or E;
    // a point closes at most one item on each of its two sides.  One item in the
    // end quadrant (q1 / q4 runs) may survive to infinity.
    let m = marks(&s.seq);
    let cnt = |a: usize, b: usize| s.seq[a..b].iter().filter(|&&x| is_item(x)).count();
    let len = s.seq.len();
    let sp1 = if endq == "q1" { spare } else { 0 };
    let sp4 = if endq == "q4" { spare } else { 0 };
    let need = match s.stage {
        // [Q2] W [Q3]: W points close one of each
        1 => cnt(0, m[0]).max(cnt(m[0] + 1, len)),
        // [Q1] N [Q2] W [Q3]: Q2 needs N, Q3 needs S, Q1 needs N or E
        2 => {
            let (r, u, l) = (cnt(0, m[0]), cnt(m[0] + 1, m[1]), cnt(m[1] + 1, len));
            l + u.max(r.saturating_sub(sp1))
        }
        // [Q1] E [Q4] S [Q3]: Q1 needs E, Q3 needs S, Q4 needs S or E
        3 => {
            let (r, q4, q3) = (cnt(0, m[0]), cnt(m[0] + 1, m[1]), cnt(m[1] + 1, len));
            (q3 + r.saturating_sub(sp1)).max(q4.saturating_sub(sp4))
        }
        _ => cnt(0, m[0]).saturating_sub(sp1).max(cnt(m[0] + 1, len).saturating_sub(sp4)),
    };
    need <= rem
}

/// Windows 11 keeps background console jobs on the E-cores (EcoQoS); opt out of that throttling, and run at
/// below-normal priority so the machine stays responsive.
#[cfg(windows)]
fn be_a_good_citizen() {
    #[repr(C)]
    struct ProcessPowerThrottlingState {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    unsafe extern "system" {
        fn GetCurrentProcess() -> isize;
        fn SetProcessInformation(h: isize, class: i32, info: *const std::ffi::c_void, size: u32) -> i32;
        fn SetPriorityClass(h: isize, class: u32) -> i32;
    }
    const PROCESS_POWER_THROTTLING: i32 = 4;
    const EXECUTION_SPEED: u32 = 0x1;
    const BELOW_NORMAL_PRIORITY_CLASS: u32 = 0x4000;
    let st = ProcessPowerThrottlingState { version: 1, control_mask: EXECUTION_SPEED, state_mask: 0 };
    unsafe {
        if SetProcessInformation(
            GetCurrentProcess(),
            PROCESS_POWER_THROTTLING,
            &st as *const _ as *const std::ffi::c_void,
            std::mem::size_of::<ProcessPowerThrottlingState>() as u32,
        ) == 0
        {
            eprintln!("note: could not disable Windows power throttling (EcoQoS)");
        }
        SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
    }
}

#[cfg(not(windows))]
fn be_a_good_citizen() {}

fn main() {
    be_a_good_citizen();
    let args: Vec<String> = std::env::args().collect();
    let nmax: usize = args[1].parse().unwrap();
    let endq = args[2].as_str();
    let mut startq = "q3";
    let mut endfirst = false;
    for a in &args[3..] {
        let forbid = |st: u8| FORBID.fetch_or(1 << st, std::sync::atomic::Ordering::Relaxed);
        match a.as_str() {
            "first" => TRACK_FIRST.store(true, std::sync::atomic::Ordering::Relaxed),
            "noW" => _ = forbid(1),
            "noN" => _ = forbid(2),
            "noS" => _ = forbid(3),
            "noE" => _ = forbid(4),
            // two ends in one quadrant: only the order with END on the cut above START
            // (START nearer the S ray in Q3, nearer the W ray in Q2)
            "endfirst" => endfirst = true,
            "start=q2" => startq = "q2",
            "start=q3" => startq = "q3",
            _ => panic!("unknown option {a}"),
        }
    }
    assert!(["q1", "q2", "q3", "q4"].contains(&endq), "END must be q1..q4");
    // Ends at infinity in Q2 / Q3 sit at x = -inf above / below the W ray from the
    // start; two ends in the same quadrant can come in either order along the rim.
    let st = |seq: Vec<u8>| State { stage: 1, done: false, first: 0, seq };
    let mut q2 = vec![START];
    let mut q3 = vec![];
    if startq == "q3" {
        q2.clear();
        q3.push(START);
    }
    let mut inits = vec![];
    let orders: Vec<(Vec<u8>, Vec<u8>)> = match endq {
        "q2" => vec![([&q2[..], &[END]].concat(), q3.clone()), ([&[END], &q2[..]].concat(), q3.clone())],
        "q3" => vec![(q2.clone(), [&q3[..], &[END]].concat()), (q2.clone(), [&[END], &q3[..]].concat())],
        _ => vec![(q2.clone(), q3.clone())],
    };
    for (a, b) in orders {
        let v = [&a[..], &[MARK], &b[..]].concat();
        let (ps, pe) = (v.iter().position(|&x| x == START), v.iter().position(|&x| x == END));
        if endfirst && matches!((ps, pe), (Some(i), Some(j)) if i < j) {
            continue;
        }
        if !inits.iter().any(|s: &State| s.seq == v) {
            inits.push(st(v));
        }
    }
    let spare = if endq == "q1" || endq == "q4" { 1 } else { 0 };
    let mut layer = Map::default();
    for s in &inits {
        *layer.entry(encode(s)).or_insert(0) += 1;
    }
    let t0 = std::time::Instant::now();
    let mut buf = Vec::new();
    for t in 0..=nmax {
        // closure under free transitions (stages only increase)
        let mut todo: Vec<(u128, u128)> = layer.iter().map(|(&k, &c)| (k, c)).collect();
        while !todo.is_empty() {
            let mut nxt = Vec::new();
            for (k, c) in todo {
                if let Some(s2) = transition(&decode(k)) {
                    let k2 = encode(&s2);
                    *layer.entry(k2).or_insert(0) += c;
                    nxt.push((k2, c));
                }
            }
            todo = nxt;
        }
        // read off finals
        let (mut tot, mut xfirst, mut xx) = (0u128, 0u128, 0u128);
        for (&k, &c) in &layer {
            let s = decode(k);
            if s.stage != 4 {
                continue;
            }
            let ok = match endq {
                "q1" => !s.done && s.seq == [START, MARK],
                "q4" => !s.done && s.seq == [MARK, START],
                _ => s.done && s.seq == [MARK],
            };
            if ok {
                tot += c;
                if s.first & 3 == 1 {
                    xfirst += c;
                }
                if s.first == 1 | 1 << 2 {
                    xx += c;
                }
            }
        }
        let mut hist = [0usize; 5];
        let mut maxlen = 0;
        for &k in layer.keys() {
            hist[(k & 7) as usize] += 1;
            maxlen = maxlen.max(decode(k).seq.len());
        }
        hist[0] = maxlen;
        println!("{} {} {} {} states {} {:.1}s stages {:?}", t, tot, xfirst, xx, layer.len(), t0.elapsed().as_secs_f64(), &hist);
        if t == nmax {
            break;
        }
        let rem = nmax - t - 1;
        let mut new = Map::default();
        for (&k, &c) in &layer {
            buf.clear();
            step(&decode(k), &mut buf);
            for s2 in &buf {
                if feasible(s2, rem, spare, endq) {
                    *new.entry(encode(s2)).or_insert(0) += c;
                }
            }
        }
        layer = new;
    }
}
