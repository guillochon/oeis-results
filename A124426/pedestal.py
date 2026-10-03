"""W_n(m): configurations when the table is replaced by m pedestals (each one a flat head).
Claim: a(n) = e^{-1} * sum_m W_n(m)/m!  (Dobinski-type)."""
from collections import defaultdict
import sympy as sp
def W(N, m):
    st = {(m, 0, 0, 0): 1}; out = []
    def place(states, kind):
        ns = defaultdict(int)
        for (s, f, e, t), w in states.items():
            opts = [(s, -1, 0, 0, 0, False)]
            opts += [(f, 0, -1, 0, 0, False), (e, 0, 0, -1, 0, True)] if kind in 'DB' else [(t, 0, 0, 0, -1, True)]
            for mu, ds, df, de, dt, enc in opts:
                if not mu: continue
                s2, f2, e2, t2 = s + ds, f + df, e + de, t + dt
                if kind == 'D': e2 += 1; s2 += (not enc)
                if kind == 'B':
                    if enc: e2 += 1
                    else: f2 += 1
                if kind == 'T': t2 += 1; s2 += (not enc)
                ns[(s2, f2, e2, t2)] += w * mu
        return ns
    for _ in range(N):
        a = place(st, 'D'); b = place(place(st, 'B'), 'T')
        st = defaultdict(int)
        for d in (a, b):
            for k, v in d.items(): st[k] += v
        out.append(sum(st.values()))
    return out
if __name__ == '__main__':
    from stack_fast import seq
    M = sp.symbols('m')
    N = 7
    tab = [W(N, m) for m in range(0, 3 * N + 3)]
    for n in range(1, N + 1):
        pts = [(m, tab[m][n - 1]) for m in range(len(tab))]
        P = sp.expand(sp.interpolate(pts, M))
        print(n, sp.factor(P), '| deg', sp.degree(P, M))
    # check Dobinski numerically
    target = seq(N)
    for n in range(1, N + 1):
        tot = sp.Rational(0)
        mm = 0; terms = []
        val = sum(sp.Rational(W(n, mm)[n - 1], sp.factorial(mm)) for mm in range(60))
        print(n, target[n - 1], sp.N(val / sp.E, 30))
