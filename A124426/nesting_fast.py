"""Fast count for nesting only (no stacking): A124426(n) = Bell(n)*Bell(n+1).
Insert dolls largest first. Closed dolls and bottom halves form chains (c = number of chains);
loose top halves form their own chains (t). Each piece starts a new chain or goes at the end of one."""
from collections import Counter
def seq(N):
    st = Counter({(0, 0): 1}); out = []
    for _ in range(N):
        new = Counter()
        for (c, t), w in st.items():
            new[c + 1, t] += w; new[c, t] += w * c              # closed doll
            for c2, m in ((c + 1, 1), (c, c)):                  # open doll: bottom, then top
                if m: new[c2, t + 1] += w * m; new[c2, t] += w * m * t
        st = new; out.append(sum(st.values()))
    return out
def bell(n):
    row = [1]
    for _ in range(n): 
        nr = [row[-1]]
        for x in row: nr.append(nr[-1] + x)
        row = nr
    return row[0]
if __name__ == '__main__':
    s = seq(20)
    assert s == [bell(n) * bell(n + 1) for n in range(1, 21)]
    print(s)
