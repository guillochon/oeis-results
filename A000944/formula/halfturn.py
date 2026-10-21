"""Half-turn-symmetric c-nets from half-turn-symmetric nonseparable maps.

All series are in M-variables: x = vertices, y = faces (networks: internal vertices/faces).
K[t][v,f] = sum over c-nets X with (v,f) of #(half-turns of pole type t)/|Aut X|, t in
VV VF FF EE EV EF (unknown, solved for).  W^B[t] = same for nonseparable maps (Liskovets-Walsh).

Tutte decomposition + dissymmetry on the subdivided tree (fixed set of a half-turn g is a path of
fixed nodes joined by side-preserving tree edges, or a single side-swapping tree edge):
  W^B = digon + sum_{fixed nodes} - sum_{side-preserving tree edges} + sum_{side-swapping tree edges}
Symmetric networks (invariant under the half-turn rho swapping the poles and the two sides) are
refined by the type of rho's fixed point inside the network ("center"): V, E or F.
"""
from fractions import Fraction

from rooted import networks
from rotations import C, Dq, E, F
from series import BS

TYPES = ('VV', 'VF', 'FF', 'EE', 'EV', 'EF')
SIG = ('V', 'E', 'F')


def tkey(a, b):
    """Canonical type name for an unordered pole pair."""
    order = {'E': 0, 'V': 1, 'F': 2}
    a, b = sorted((a, b), key=lambda c: order[c])
    s = a + b
    return {'VE': 'EV', 'FE': 'EF', 'FV': 'VF'}.get(s, s)


def WB(N):
    """Liskovets-Walsh weights for half-turns of nonseparable maps, by pole type."""
    out = {t: BS(N) for t in TYPES}
    for v in range(2, N + 1):
        for f in range(2, N + 1 - v):
            i, j = v - 1, f - 1
            e = i + j
            q = {}
            if i % 2 == 0 and j % 2 == 0:
                q['VF'] = E(i // 2, j // 2)
            if i % 2 == 1:
                q['VV'] = F((i - 1) // 2, (j + 1) // 2)
            if j % 2 == 1:
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


def geom(N, a):
    """a/(1-a)."""
    return a * (1 - a).inv()


def solve(N, iters=None, verbose=False):
    D, S, P, H = networks(N)
    x = BS.mono(N, 1, 0)
    y = BS.mono(N, 0, 1)
    A = D - S
    Bn = D - P
    D2, A2, Bn2, H2 = D.scale(2, 2), A.scale(2, 2), Bn.scale(2, 2), H.scale(2, 2)
    D2pow = [BS.const(N, 1)]
    for _ in range(N):
        D2pow.append(D2pow[-1] * D2)
    wb = WB(N)
    xy = x * y
    digon = {'VV': xy * xy * Fraction(1, 8), 'FF': xy * xy * Fraction(1, 8), 'EE': xy * xy * Fraction(1, 8)}

    # node/edge terms that do not involve K except through symmetric networks
    gA = geom(N, x * x * A2)          # sum_{m>=1} (x^2 A2)^m
    gB = geom(N, y * y * Bn2)
    K = {t: BS(N) for t in TYPES}

    def base(t, k):
        """sum_X K x^v y^f D2^(free edge pairs), with root-adjusted exponents k = (dv, df, pe)."""
        dv, df, pe = k
        out = BS(N)
        for v, f, c in K[t].items():
            e = v + f - 2
            fp = (e - pe) // 2
            if v + dv < 0 or f + df < 0:
                continue
            out = out + (D2pow[fp] * c).shift(v + dv, f + df) if fp <= N else out
        return out

    for it in range(iters or (N + 2)):
        # symmetric networks, fixed point
        Hrev_EV = base('EV', (-2, -2, 1)) * 4
        Hrev_EF = base('EF', (-2, -2, 1)) * 4
        Hrev_EE = base('EE', (-2, -2, 2)) * 8
        Srev = {s: BS(N) for s in SIG}
        Prev = {s: BS(N) for s in SIG}
        Hrev = {s: BS(N) for s in SIG}
        for _ in range(N + 2):
            Drev = {s: (1 if s == 'E' else 0) + Srev[s] + Prev[s] + Hrev[s] for s in SIG}
            Arev = {s: (1 if s == 'E' else 0) + Prev[s] + Hrev[s] for s in SIG}
            Bnrev = {s: (1 if s == 'E' else 0) + Srev[s] + Hrev[s] for s in SIG}
            nS = {s: gA * Arev[s] for s in SIG}
            nS['V'] = nS['V'] + x * A2 * (1 - x * x * A2).inv()
            nP = {s: gB * Bnrev[s] for s in SIG}
            nP['F'] = nP['F'] + y * Bn2 * (1 - y * y * Bn2).inv()
            nH = {s: Hrev_EE * Drev[s] for s in SIG}
            nH['V'] = nH['V'] + Hrev_EV
            nH['F'] = nH['F'] + Hrev_EF
            if all(nS[s].c == Srev[s].c and nP[s].c == Prev[s].c and nH[s].c == Hrev[s].c for s in SIG):
                break
            Srev, Prev, Hrev = nS, nP, nH
        Drev = {s: (1 if s == 'E' else 0) + Srev[s] + Prev[s] + Hrev[s] for s in SIG}
        Arev = {s: (1 if s == 'E' else 0) + Prev[s] + Hrev[s] for s in SIG}
        Bnrev = {s: (1 if s == 'E' else 0) + Srev[s] + Hrev[s] for s in SIG}

        rhs = {t: BS(N) for t in TYPES}
        for t in ('VV', 'VF', 'FF'):
            rhs[t] = rhs[t] + base(t, (0, 0, 0))
        bEV = base('EV', (0, 0, 1))
        bEF = base('EF', (0, 0, 1))
        bEE = base('EE', (0, 0, 2))
        for s in SIG:
            rhs[tkey(s, 'V')] = rhs[tkey(s, 'V')] + bEV * Drev[s]
            rhs[tkey(s, 'F')] = rhs[tkey(s, 'F')] + bEF * Drev[s]
            for s2 in SIG:
                rhs[tkey(s, s2)] = rhs[tkey(s, s2)] + bEE * Drev[s] * Drev[s2]
        # S-node (polygon, l >= 3)
        for m in range(2, N + 1):
            term = (x * x * A2)
            pw = BS.const(N, 1)
            for _ in range(m):
                pw = pw * term
            if not any(any(r) for r in pw.c):
                break
            rhs['FF'] = rhs['FF'] + (y * y * pw) * Fraction(1, 8 * m)
            rhs['VV'] = rhs['VV'] + (y * y * pw) * Fraction(1, 8)
        # edge-edge (l = 2m, m >= 2): (1/8) x^2m y^2 A2^(m-1) Arev Arev = (1/8) y^2 x^2 (x^2A2)^(m-1) ...
        ee = y * y * x * x * gA * Fraction(1, 8)    # sum_{m>=2} x^{2m} A2^{m-1} = x^2 * sum_{m>=1} (x^2 A2)^m
        for s in SIG:
            for s2 in SIG:
                rhs[tkey(s, s2)] = rhs[tkey(s, s2)] + ee * Arev[s] * Arev[s2]
            # vertex-edge (l = 2m+1, m >= 1): (1/4) x^(2m+1) y^2 A2^m Arev
            rhs[tkey('V', s)] = rhs[tkey('V', s)] + x * y * y * gA * Arev[s] * Fraction(1, 4)
        # P-node (bond, l >= 3), dual
        for m in range(2, N + 1):
            pw = BS.const(N, 1)
            for _ in range(m):
                pw = pw * (y * y * Bn2)
            if not any(any(r) for r in pw.c):
                break
            rhs['VV'] = rhs['VV'] + (x * x * pw) * Fraction(1, 8 * m)
            rhs['FF'] = rhs['FF'] + (x * x * pw) * Fraction(1, 8)
        ee = x * x * y * y * gB * Fraction(1, 8)
        for s in SIG:
            for s2 in SIG:
                rhs[tkey(s, s2)] = rhs[tkey(s, s2)] + ee * Bnrev[s] * Bnrev[s2]
            rhs[tkey('F', s)] = rhs[tkey('F', s)] + x * x * y * gB * Bnrev[s] * Fraction(1, 4)
        # side-preserving tree edges (subtract): digon skeleton, one EE half-turn, |Aut| = 8
        net = {'S': Srev, 'P': Prev, 'H': Hrev}
        for t1 in 'SPH':
            for t2 in 'SPH':
                if t1 == t2 and t1 in 'SP':
                    continue
                for s in SIG:
                    for s2 in SIG:
                        rhs[tkey(s, s2)] = rhs[tkey(s, s2)] - xy * net[t1][s] * net[t2][s2] * Fraction(1, 8)
        # side-swapping tree edges (add): H-networks swapped, poles VV or FF
        rhs['VV'] = rhs['VV'] + xy * H2 * Fraction(1, 8)
        rhs['FF'] = rhs['FF'] + xy * H2 * Fraction(1, 8)
        for t, d in digon.items():
            rhs[t] = rhs[t] + d
        # update
        done = True
        for t in TYPES:
            res = wb[t] - rhs[t]
            for v, f, c in res.items():
                done = False
                K[t][v, f] = K[t][v, f] + c
        if verbose:
            print('iter', it, 'done' if done else '')
        if done:
            break
    return K


def main():
    import json
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    K = solve(N, verbose=True)
    stats = {}
    for fn in sys.argv[2:]:
        for key, v in json.load(open(fn)).items():
            stats.setdefault(key, {}).update({t: Fraction(w) for t, w in v.items()})
    bad = checked = 0
    for key, d in stats.items():
        v, f = map(int, key.split(','))
        if v + f > N:
            continue
        for t in TYPES:
            pol = ''.join(sorted(t))
            want = d.get(json.dumps(['rot', 2, pol]), 0)
            got = K[t][v, f]
            checked += 1
            if got != want:
                bad += 1
                print('MISMATCH', t, v, f, got, want)
    print('checked', checked, 'bad', bad)


if __name__ == '__main__':
    main()
