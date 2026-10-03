from collections import Counter
def a(n):
    # state (s,f,e,t): free heads, empty open bottoms in towers, outer closed dolls, outer tops
    st = Counter({(0, 0, 0, 0): 1})
    for _ in range(n):
        new = Counter()
        for (s, f, e, t), w in st.items():
            # closed doll: table, on a head, in an open bottom, inside a closed doll
            new[s+1, f, e+1, t] += w; new[s, f, e+1, t] += w*s
            new[s+1, f-1, e+1, t] += w*f; new[s, f, e, t] += w*e
            # open doll: bottom (table, head, open bottom, closed doll), then top (table, head, top)
            for (s2, f2), m in (((s, f+1), 1), ((s-1, f+1), s), ((s, f), f+e)):
                if m:
                    new[s2+1, f2, e, t+1] += w*m; new[s2, f2, e, t+1] += w*m*s2
                    new[s2, f2, e, t] += w*m*t
        st = new
    return sum(st.values())
print([a(n) for n in range(21)])
