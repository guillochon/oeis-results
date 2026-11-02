"""Transfer-matrix count of meanders crossing two perpendicular roads (prototype).

The processed region is the part left of a cut curve that crosses some of the rays
W, N, E, S of the axes.  The cut is a sequence (top to bottom) of markers (where it
crosses a ray) and items (river edges crossing the cut, each labelled by what its
processed end connects to: START, END, or the partner item of a pair).
Stages, each processing the points of one ray in a fixed direction:
  A  cut x = c < 0, marker W moves from -inf towards O         seq: U  W  L
  B  marker N inserted at the top, moves down from +inf to O    seq: R  N  U' W  L
  -- N and W adjacent at O: push the cut across O:  N W -> E S
  C  marker S moves from O down to -inf                         seq: R  E  Q4 S  Q3
  -- nothing after S: drop it
  D  cut x = c > 0, marker E moves from O to +inf               seq: Q1 E  Q4
Processing a point at the active marker: each of its two half-edges either closes
the item adjacent to the marker on its side or opens a new adjacent item
(Jensen's meander step).  Closed loops are forbidden.
"""
import sys
from collections import defaultdict

W, N, E, S = "W", "N", "E", "S"
START, END = "s", "e"
ACTIVE = {1: W, 2: N, 3: S, 4: E}


def norm(seq):
    """Relabel pair ids by first appearance."""
    m, out, k = {}, [], 0
    for x in seq:
        if isinstance(x, int):
            if x not in m:
                m[x] = k
                k += 1
            out.append(m[x])
        else:
            out.append(x)
    return tuple(out)


def stage(seq):
    ms = tuple(x for x in seq if x in (W, N, E, S))
    return {(W,): 1, (N, W): 2, (E, S): 3, (E,): 4}[ms]


def is_item(x):
    return x not in (W, N, E, S)


def step(seq, done):
    """All states after processing one point at the active marker."""
    if done:
        return []
    st = stage(seq)
    i = seq.index(ACTIVE[st])
    out = []
    before = i > 0 and is_item(seq[i - 1])
    after = i + 1 < len(seq) and is_item(seq[i + 1])
    fresh = max([x for x in seq if isinstance(x, int)], default=-1) + 1
    # open both: new pair around the marker
    out.append((seq[:i] + (fresh, seq[i], fresh) + seq[i + 1:], False))
    # close before, open after: edge continues
    if before:
        out.append((seq[:i - 1] + (seq[i], seq[i - 1]) + seq[i + 1:], False))
    if after:
        out.append((seq[:i] + (seq[i + 1], seq[i]) + seq[i + 2:], False))
    if before and after:
        a, b = seq[i - 1], seq[i + 1]
        rest = list(seq[:i - 1] + (seq[i],) + seq[i + 2:])
        if a == b and isinstance(a, int):
            pass                                   # closed loop
        elif {a, b} == {START, END}:
            if not any(is_item(x) for x in rest):
                out.append((tuple(rest), True))
        elif isinstance(a, int) and isinstance(b, int):
            j = [k for k, x in enumerate(rest) if x in (a, b)]
            # the partners of a and b become a pair
            for k in j:
                rest[k] = a
            out.append((tuple(rest), False))
        else:                                      # one of them is START/END
            lab, other = (a, b) if not isinstance(a, int) else (b, a)
            if isinstance(other, int):
                rest = [lab if x == other else x for x in rest]
                out.append((tuple(rest), False))
            # START-START or END-END impossible
    return [(norm(s), d) for s, d in out]


def transitions(seq, done):
    """Free stage changes (no points)."""
    st = stage(seq)
    if st == 1:
        return [((N,) + seq, done)]
    if st == 2:
        i = seq.index(N)
        if seq[i + 1] == W:
            return [(seq[:i] + (E, S) + seq[i + 2:], done)]
    if st == 3:
        i = seq.index(S)
        if i == len(seq) - 1:
            return [(seq[:i], done)]
    return []


def count(n, end_quad):
    """end_quad: 'Q3', 'Q4', 'Q1', 'Q2'.  Returns the number of meanders with n crossings."""
    if end_quad == "Q3":
        inits = [(W, START, END), (W, END, START)]
    elif end_quad == "Q2":
        inits = [(END, W, START)]
    else:
        inits = [(W, START)]
    layer = defaultdict(int)
    for s in inits:
        layer[(norm(s), False)] += 1
    for t in range(n + 1):
        # closure under free transitions
        todo = list(layer.items())
        while todo:
            nxt = []
            for (s, d), c in todo:
                for s2, d2 in transitions(s, d):
                    layer[(s2, d2)] += c
                    nxt.append(((s2, d2), c))
            todo = nxt
        if t == n:
            break
        new = defaultdict(int)
        for (s, d), c in layer.items():
            items = sum(1 for x in s if is_item(x))
            if items > 2 * (n - t) + 1:            # cannot all be closed in time
                continue
            for s2, d2 in step(s, d):
                new[(s2, d2)] += c
        layer = new
    tot = 0
    for (s, d), c in layer.items():
        if stage(s) != 4:
            continue
        if end_quad == "Q1" and s == (START, E) and not d:
            tot += c
        elif end_quad == "Q4" and s == (E, START) and not d:
            tot += c
        elif end_quad in ("Q2", "Q3") and s == (E,) and d:
            tot += c
    return tot


if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    for n in range(M + 1):
        per = [count(n, q) for q in ("Q3", "Q4", "Q1", "Q2")]
        a906 = per[2] if n % 2 == 0 else per[1]
        print(n, a906, sum(per), per, flush=True)
