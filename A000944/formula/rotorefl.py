"""Rotatory reflections without a reflection in <g>: first the antipodal map (g^2 = 1).

All series in half variables U = x^2, W = y^2 (an antipodally symmetric map has v = 2a, f = 2b).
1. General planar maps M with the antipodal map <-> maps on the projective plane (the quotient):
   W^gen(a, b) = Qrp2(a, b) / (4e), e = 2a + 2b - 2 (flag-rooted pairs (M, g) = rooted quotients).
2. Blocks: g fixes a unique block B (no fixed vertex, so no fixed cut vertex; a bridge or a loop would
   contain a fixed point). Its 2 e_B corners form e_B free orbits, each carrying a sequence of pendant
   rooted map Gamma = P(1)/(xy):  W^gen = sum K^B(a,b) U^a W^b (Gamma^2)^(a+b-1).
3. Tutte decomposition: g fixes a unique node or inverts one tree edge (then both sides are swapped
   h-networks):  W^B = sum K^C U^a W^b D^(a+b-1) + polygons + bonds + digon + inverted edge.
"""
from fractions import Fraction

from ms import Ctx
from pipeline import invert_D, networks, powers
from rp2 import maps_sphere_rp2
from t3 import T3


def antipodal(ctx, P=None, Q=None, T=None):
    """K^C(a, b): c-net weights for the antipodal map, half variables, as an MS series."""
    N, p = ctx.N, ctx.p
    if T is None:
        T = T3(N + 2, N + 2, 2 * N + 6, p)
        P, Q = maps_sphere_rp2(T)
    G = ctx.zero()
    for a in range(1, N + 1):
        for b in range(1, N + 1 - a):
            if a > T.X or b > T.Y:
                continue
            q = int(Q[a, b].sum() % p)
            if q:
                e = 2 * a + 2 * b - 2
                G[a, b] = Fraction(q, 4 * e) if False else q * pow(4 * e, p - 2, p) % p
    # pendant sequences at a corner, in half variables
    R1 = ctx.zero()
    for a in range(N + 1):
        for b in range(N + 1 - a):
            if a + 1 <= T.X and b + 1 <= T.Y:
                R1[a, b] = int(P[a + 1, b + 1].sum() % p)
    Gam = R1                             # everything hanging in one corner: a rooted map, P(1) / (x y)
    G2 = Gam * Gam
    KB = invert_D(G * G2, powers(G2, N + 1))     # sum K U^a W^b G2^(a+b) = G * G2
    D, S, Pn, H = networks(ctx)
    U, Wv = ctx.mono(1, 0), ctx.mono(0, 1)
    A, Bn = D - S, D - Pn
    rest = ctx.zero()
    pa = ctx.const(1)
    pb = ctx.const(1)
    for m in range(1, N + 1):
        pa = pa * (U * A)
        pb = pb * (Wv * Bn)
        if m >= 2:
            rest = rest + (Wv * pa + U * pb) * Fraction(1, 8 * m)
    rest = rest + U * Wv * (H + 1) * Fraction(1, 8)      # inverted tree edge + digon
    KC = invert_D((KB - rest) * D, powers(D, N + 1))      # sum K^C U^a W^b D^(a+b) = (...) * D
    return KC


if __name__ == '__main__':
    import glob
    import json
    import sys
    from ms import PRIMES, frac_mod
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    p = PRIMES[0]
    ctx = Ctx(N, p)
    KC = antipodal(ctx)
    stats = {}
    for fn in glob.glob('../brute/stats_*.json'):
        for key, d in json.load(open(fn)).items():
            stats.setdefault(key, {}).update(d)
    ok = bad = 0
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v + f > 2 * N:
            continue
        want = Fraction(d.get('["rr", 2]', '0'))
        got = KC[v // 2, f // 2] if v % 2 == 0 and f % 2 == 0 else 0
        if got == frac_mod(want, p):
            ok += 1
        else:
            bad += 1
            print('mismatch', v, f, 'brute', want)
    print('ok', ok, 'bad', bad)


# ------------------------------------------------------------------ order n >= 4 (n even)
#
# g of even order n >= 4 with no reflection in <g>: g^2 is a rotation of order n/2 whose poles p, q are
# swapped by g (both vertices, both faces, or - for n = 4 - both edge midpoints).
# 1. Quotient: projective plane with one cone point of order n/2 at a vertex, a face, or (n = 4) the end
#    of a semi-edge. pi_1 of the orbifold is cyclic of order n, so each rooted quotient with its cone
#    cell lifts to phi(n) pairs (M, g):
#      VV: v = n(vQ - 1) + 2, f = n fQ:  W = phi(n) R_Q vQ / (4e);   FF symmetric;
#      EE (n = 4): v = 4 vQ, f = 4 fQ, semi-edge in one of 2 eQ corners:  W = R_Q / 4.
# 2. Blocks: the fixed block's 2e corners form free orbits of size n: Gamma(x^n, y^n)^(2e/n).
# 3. Tutte decomposition: unique fixed node. Edge orbits have size n, except for n = 4 with edge
#    poles: the two pole edges (orbit of size 2) carry half-turn-symmetric networks (Drev, by centre).


def _mvar(ctx, quot, n, dv, df):
    """quotient series (U = x^n, W = y^n) -> M variables, times x^dv y^df."""
    return quot.scale(n, n).shift(dv, df)


def even_order(ctx, n, nets, Drev, P, Q, T):
    """{'VV','FF'(,'EE')} -> M-variable series of c-net weights for rotatory reflections of order n."""
    from math import gcd
    N, p = ctx.N, ctx.p
    phin = sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)
    inv = lambda k: pow(k % p, p - 2, p)
    types = ('VV', 'FF', 'EE') if n == 4 else ('VV', 'FF')
    Rq = {}
    for a in range(1, T.X + 1):
        for b in range(1, T.Y + 1):
            r = int(Q[a, b].sum() % p)
            if r:
                Rq[(a, b)] = r
    # 1. general maps, M variables
    Wg = {t: ctx.zero() for t in types}
    for (vq, fq), r in Rq.items():
        eq = vq + fq - 1
        e = n * eq
        for t in types:
            if t == 'VV':
                v, f, w = n * (vq - 1) + 2, n * fq, phin * r * vq % p * inv(4 * e) % p
            elif t == 'FF':
                v, f, w = n * vq, n * (fq - 1) + 2, phin * r * fq % p * inv(4 * e) % p
            else:
                v, f, w = 4 * vq, 4 * fq, r * inv(4) % p
            if v + f <= N:
                Wg[t].a[v, f] = (Wg[t].a[v, f] + w) % p
    # corner decoration Gamma(x^n, y^n) and network series at (x^n, y^n)
    G1 = ctx.zero()
    for a in range(N + 1):
        for b in range(N + 1 - a):
            if a + 1 <= T.X and b + 1 <= T.Y:
                G1[a, b] = int(P[a + 1, b + 1].sum() % p)
    Gn = G1.scale(n, n)
    D, S, Pn, H = nets
    Dn, An, Bnn = D.scale(n, n), (D - S).scale(n, n), (D - Pn).scale(n, n)
    x, y = ctx.mono(1, 0), ctx.mono(0, 1)
    Gpow = powers(Gn, 2 * N // n + 2)
    Dpow = powers(Dn, 2 * N // n + 2)

    def apply(K, pw, expo):
        out = ctx.zero()
        for v, f, c in K.items():
            k = expo(v, f)
            if 0 <= k < len(pw):
                out = out + pw[k].shift(v, f) * c
        return out

    corners = lambda v, f: 2 * (v + f - 2) // n
    # 2. invert the block decomposition (triangular: leading coefficient 1)
    KB = {t: ctx.zero() for t in types}
    for _ in range(N + 2):
        done = True
        for t in types:
            res = Wg[t] - apply(KB[t], Gpow, corners)
            if not res.iszero():
                done = False
                KB[t] = KB[t] + res
        if done:
            break
    # 3. invert the Tutte decomposition
    fixed = {t: ctx.zero() for t in types}
    pa, pb = ctx.const(1), ctx.const(1)
    for k in range(1, N + 1):
        pa, pb = pa * (x.scale(n, 1) * An), pb * (y.scale(1, n) * Bnn)   # (x^n A_n)^k, (y^n B_n)^k
        if pa.iszero() and pb.iszero():
            break
        fixed['FF'] = fixed['FF'] + (y * y * pa) * (phin * inv(4 * n * k) % p)
        fixed['VV'] = fixed['VV'] + (x * x * pb) * (phin * inv(4 * n * k) % p)
    free = lambda v, f: (v + f - 2) // n                   # VV, FF: all edges in orbits of size n
    freeE = lambda v, f: (v + f - 4) // n                  # EE: two pole edges in an orbit of size 2
    sig2 = {'V': 'VV', 'E': 'EE', 'F': 'FF'}
    KC = {t: ctx.zero() for t in types}
    for _ in range(N + 2):
        rhs = {t: fixed[t] for t in types}
        rhs['VV'] = rhs['VV'] + apply(KC['VV'], Dpow, free)
        rhs['FF'] = rhs['FF'] + apply(KC['FF'], Dpow, free)
        if n == 4:
            base = apply(KC['EE'], Dpow, freeE)
            for s_, t in sig2.items():
                rhs[t] = rhs[t] + base * Drev[s_].scale(2, 2)
        done = True
        for t in types:
            res = KB[t] - rhs[t]
            if not res.iszero():
                done = False
                KC[t] = KC[t] + res
        if done:
            return KC
    raise RuntimeError('even-order rotatory reflections did not converge')
