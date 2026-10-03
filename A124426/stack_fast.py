"""Fast count: insert dolls largest first; each piece takes one free slot (or the table).
State: s = free flat tops, f = empty open bottoms that stick out, e = empty enclosed D/B insides,
t = empty top insides."""
from collections import defaultdict
def seq(N):
    st = {(0, 0, 0, 0): 1}; out = []
    for _ in range(N):
        def place(states, kind):
            ns = defaultdict(int)
            for (s, f, e, t), w in states.items():
                opts = []  # (multiplicity, ds, df, de, dt, enclosed)
                if kind in 'DB':
                    opts += [(1, 0, 0, 0, 0, False), (s, -1, 0, 0, 0, False), (f, 0, -1, 0, 0, False), (e, 0, 0, -1, 0, True)]
                else:
                    opts += [(1, 0, 0, 0, 0, False), (s, -1, 0, 0, 0, False), (t, 0, 0, 0, -1, True)]
                for m, ds, df, de, dt, enc in opts:
                    if not m: continue
                    s2, f2, e2, t2 = s + ds, f + df, e + de, t + dt
                    if kind == 'D': e2 += 1; s2 += (not enc)
                    if kind == 'B':
                        if enc: e2 += 1
                        else: f2 += 1
                    if kind == 'T': t2 += 1; s2 += (not enc)
                    ns[(s2, f2, e2, t2)] += w * m
            return ns
        closed = place(st, 'D')
        opened = place(place(st, 'B'), 'T')
        st = defaultdict(int)
        for d in (closed, opened):
            for k, v in d.items(): st[k] += v
        out.append(sum(st.values()))
    return out
if __name__ == '__main__':
    print(seq(20))
