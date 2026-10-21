"""Modular (Z/p) pipeline: rooted c-nets, rotations of order >= 3, half-turns.

Same equations as rooted.py / rotations.py / halfturn.py (see those docstrings), on ms.MS series.
"""
from fractions import Fraction
from math import gcd

from ms import Ctx
from rooted import B
from rotations import C, Dq, E, F


def phi(n):
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


def networks(ctx):
    N = ctx.N
    D = ctx.zero()
    for i in range(1, N + 2):
        for j in range(1, N + 3 - i):
            if (i - 1) + (j - 1) <= N:
                D[i - 1, j - 1] = B(i, j)
    x, y = ctx.mono(1, 0), ctx.mono(0, 1)
    # S = x A^2 / (1 - x A) with A = D - S gives D = A / (1 - x A), so S = x D^2 / (1 + x D); same for P.
    S = x * D * D * (1 + x * D).inv()
    P = y * D * D * (1 + y * D).inv()
    H = D - 1 - S - P
    return D, S, P, H


def powers(f, n):
    out = [f.ctx.const(1)]
    for _ in range(n):
        out.append(out[-1] * f)
    return out


def invert_D(G, Dpow):
    """m with sum_{a,b} m_ab u^a w^b D^(a+b) = G (u, w the series variables)."""
    ctx = G.ctx
    N = ctx.N
    m = ctx.zero()
    R = G.a.copy()
    p = ctx.p
    for d in range(N + 1):
        for i in range(d + 1):
            j = d - i
            r = int(R[i, j])
            if r:
                m.a[i, j] = r
                R = (R - r * Dpow[d].shift(i, j).a) % p
    return m


def subst_D(m, Dpow):
    """sum_{a,b} m_ab u^a w^b D^(a+b)."""
    ctx = m.ctx
    out = ctx.zero()
    for a, b, c in m.items():
        out = out + Dpow[a + b].shift(a, b) * c
    return out


def rooted_cnets(ctx, nets=None):
    D, S, P, H = nets or networks(ctx)
    Dpow = powers(D, ctx.N)
    return invert_D(H * D.inv(), Dpow)


def neglog1m(A):
    ctx = A.ctx
    out = ctx.zero()
    pw = ctx.const(1)
    for n in range(1, ctx.N + 1):
        pw = pw * A
        if pw.iszero():
            break
        out = out + pw * Fraction(1, n)
    return out


def rot_big(ctx, k, tau, nets):
    """c_tau(v', f') (quotient orbit counts) for rotations of order k >= 3."""
    D, S, P, H = nets
    N = ctx.N
    L = ctx.zero()
    for a in range(N + 1):
        for b in range(N + 1 - a):
            if tau == 'VF':
                if a < 1 or b < 1:
                    continue
                q = E(a, b)
            elif tau == 'VV':
                if b < 1:
                    continue
                q = F(a, b)
            else:
                if a < 1:
                    continue
                q = F(b, a)
            L[a, b] = phi(k) * Fraction(q) / (4 * k * (a + b))
    u, w = ctx.mono(1, 0), ctx.mono(0, 1)
    if tau == 'FF':
        L = L - neglog1m(u * (1 + P + H)) * Fraction(phi(k), 4 * k)
    if tau == 'VV':
        L = L - neglog1m(w * (1 + S + H)) * Fraction(phi(k), 4 * k)
    return invert_D(L, powers(D, N))


# ---------------------------------------------------------------- half-turns

TYPES = ('VV', 'VF', 'FF', 'EE', 'EV', 'EF')
SIG = ('V', 'E', 'F')


def tkey(a, b):
    order = {'E': 0, 'V': 1, 'F': 2}
    a, b = sorted((a, b), key=lambda c: order[c])
    s = a + b
    return {'VE': 'EV', 'FE': 'EF', 'FV': 'VF'}.get(s, s)


def WB_half(ctx):
    N = ctx.N
    out = {t: ctx.zero() for t in TYPES}
    for v in range(2, N + 1):
        for f in range(2, N + 1 - v):
            i, j = v - 1, f - 1
            e = i + j
            q = {}
            if i % 2 == 0 and j % 2 == 0:
                q['VF'] = E(i // 2, j // 2)
            if i % 2 == 1 and j % 2 == 1:   # d = 2 must divide e = i + j
                q['VV'] = F((i - 1) // 2, (j + 1) // 2)
                q['FF'] = F((j - 1) // 2, (i + 1) // 2)
            if i % 2 == 1 and j % 2 == 1:
                q['EE'] = e * C(i, j)
            if i % 2 == 0 and j % 2 == 1:
                q['EV'] = e * Dq(i, j)
            if i % 2 == 1 and j % 2 == 0:
                q['EF'] = e * Dq(j, i)
            for t, val in q.items():
                out[t][v, f] = Fraction(val) / (4 * e)
    return out


def halfturn(ctx, nets, verbose=False, want_rev=False):
    """K[t]: c-net half-turn weights by pole type, M-variables (v, f)."""
    N = ctx.N
    D, S, P, H = nets
    x, y = ctx.mono(1, 0), ctx.mono(0, 1)
    xy = x * y
    A, Bn = D - S, D - P
    D2, A2, Bn2, H2 = D.scale(2, 2), A.scale(2, 2), Bn.scale(2, 2), H.scale(2, 2)
    D2pow = powers(D2, N // 2 + 1)
    gA = (x * x * A2) * (1 - x * x * A2).inv()
    gB = (y * y * Bn2) * (1 - y * y * Bn2).inv()
    wb = WB_half(ctx)
    one = ctx.const(1)
    # K-independent part of the right-hand side
    fixed = {t: ctx.zero() for t in TYPES}
    for t in ('VV', 'FF', 'EE'):
        fixed[t] = fixed[t] + xy * xy * Fraction(1, 8)            # digon
    pw = ctx.const(1)
    for m in range(1, N + 1):
        pw = pw * (x * x * A2)
        if pw.iszero():
            break
        if m >= 2:
            fixed['FF'] = fixed['FF'] + y * y * pw * Fraction(1, 8 * m)
            fixed['VV'] = fixed['VV'] + y * y * pw * Fraction(1, 8)
    pw = ctx.const(1)
    for m in range(1, N + 1):
        pw = pw * (y * y * Bn2)
        if pw.iszero():
            break
        if m >= 2:
            fixed['VV'] = fixed['VV'] + x * x * pw * Fraction(1, 8 * m)
            fixed['FF'] = fixed['FF'] + x * x * pw * Fraction(1, 8)
    fixed['VV'] = fixed['VV'] + xy * xy * H2 * Fraction(1, 8)
    fixed['FF'] = fixed['FF'] + xy * xy * H2 * Fraction(1, 8)
    sA = x * A2 * (1 - x * x * A2).inv()
    sB = y * Bn2 * (1 - y * y * Bn2).inv()

    K = {t: ctx.zero() for t in TYPES}

    def base(t, dv, df, pe):
        out = ctx.zero()
        for v, f, c in K[t].items():
            fp = (v + f - 2 - pe) // 2
            out = out + D2pow[fp].shift(v + dv, f + df) * c
        return out

    for it in range(N + 3):
        hEV = base('EV', -2, -2, 1) * 4
        hEF = base('EF', -2, -2, 1) * 4
        hEE = base('EE', -2, -2, 2) * 8
        Srev = {s: ctx.zero() for s in SIG}
        Prev = {s: ctx.zero() for s in SIG}
        Hrev = {s: ctx.zero() for s in SIG}
        for _ in range(N + 3):
            e_ = {s: (one if s == 'E' else ctx.zero()) for s in SIG}
            Drev = {s: e_[s] + Srev[s] + Prev[s] + Hrev[s] for s in SIG}
            Arev = {s: e_[s] + Prev[s] + Hrev[s] for s in SIG}
            Bnrev = {s: e_[s] + Srev[s] + Hrev[s] for s in SIG}
            nS = {s: gA * Arev[s] for s in SIG}
            nS['V'] = nS['V'] + sA
            nP = {s: gB * Bnrev[s] for s in SIG}
            nP['F'] = nP['F'] + sB
            nH = {s: hEE * Drev[s] for s in SIG}
            nH['V'] = nH['V'] + hEV
            nH['F'] = nH['F'] + hEF
            same = all(nS[s] == Srev[s] and nP[s] == Prev[s] and nH[s] == Hrev[s] for s in SIG)
            Srev, Prev, Hrev = nS, nP, nH
            if same:
                break
        e_ = {s: (one if s == 'E' else ctx.zero()) for s in SIG}
        Drev = {s: e_[s] + Srev[s] + Prev[s] + Hrev[s] for s in SIG}
        Arev = {s: e_[s] + Prev[s] + Hrev[s] for s in SIG}
        Bnrev = {s: e_[s] + Srev[s] + Hrev[s] for s in SIG}

        rhs = {t: fixed[t] for t in TYPES}
        for t in ('VV', 'VF', 'FF'):
            rhs[t] = rhs[t] + base(t, 0, 0, 0)
        bEV, bEF, bEE = base('EV', 0, 0, 1), base('EF', 0, 0, 1), base('EE', 0, 0, 2)
        eeS = x * x * y * y * gA * Fraction(1, 8)
        eeP = x * x * y * y * gB * Fraction(1, 8)
        veS = x * y * y * gA * Fraction(1, 4)
        feP = x * x * y * gB * Fraction(1, 4)
        net = {'S': Srev, 'P': Prev, 'H': Hrev}
        for s in SIG:
            rhs[tkey(s, 'V')] = rhs[tkey(s, 'V')] + bEV * Drev[s] + veS * Arev[s]
            rhs[tkey(s, 'F')] = rhs[tkey(s, 'F')] + bEF * Drev[s] + feP * Bnrev[s]
            for s2 in SIG:
                t = tkey(s, s2)
                rhs[t] = rhs[t] + bEE * Drev[s] * Drev[s2] + eeS * Arev[s] * Arev[s2] + eeP * Bnrev[s] * Bnrev[s2]
                for t1 in 'SPH':
                    for t2 in 'SPH':
                        if t1 == t2 and t1 in 'SP':
                            continue
                        rhs[t] = rhs[t] - xy * xy * net[t1][s] * net[t2][s2] * Fraction(1, 8)
        done = True
        for t in TYPES:
            res = wb[t] - rhs[t]
            if not res.iszero():
                done = False
                K[t] = K[t] + res
        if verbose:
            print('iter', it, flush=True)
        if done:
            return (K, Drev) if want_rev else K
    raise RuntimeError('half-turn iteration did not converge')
