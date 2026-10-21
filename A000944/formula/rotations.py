"""Rotation-symmetric c-nets from rotation-symmetric nonseparable maps (Liskovets-Walsh quotients).

For a rotation g of order k with pole types tau in {VF, VV, FF} (k >= 3; edge poles need k = 2),
the weighted count  W_tau(v,f) = sum_M #(elements of type tau in Aut M) / |Aut M|  (Aut includes
orientation-reversing maps) satisfies, for nonseparable maps M,
    W^B_tau = R-part + (S-part if FF) + (P-part if VV)
because for k >= 3 no tree edge of the Tutte decomposition is g-invariant (the fixed node is unique).
  R-part = x^pv y^pf * C_tau(u D(u,w), w D(u,w)),  u = x^k, w = y^k,
  C_tau(X,Y) = sum c_tau(v',f') X^v' Y^f'  with c_tau the c-net weights (v = k v' + pv, f = k f' + pf),
  S-part = phi(k)/(4k) * -log(1 - u (1+P+H)(u,w))   [polygon, poles in its two faces],
  P-part = phi(k)/(4k) * -log(1 - w (1+S+H)(u,w))   [bond, poles at its two vertices].
W^B_tau = phi(k) Q_tau / (4e) with Q_tau the dart-rooted quotient counts of Liskovets-Walsh
(A342061 program of A. Howroyd: functions E and F below).
"""
from fractions import Fraction
from math import comb, gcd

from rooted import B, invert_xy, networks
from series import BS


def phi(n):
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


def C(i, j):
    return Fraction((i + j - 1) * (i + 1) * (j + 1), 2 * (2 * i + j - 1) * (2 * j + i - 1)) * B((i + 1) // 2, (j + 1) // 2)


def Dq(i, j):
    return Fraction(j + 1, 2) * B(i // 2, (j + 1) // 2)


def E(i, j):
    return ((i - 1) * (j - 1) + 2 * (i + j) * (i + j - 1)) * B(i, j)


def F(i, j):
    # PARI: if(!i, j==1, main*B[i,j]) + if(j-1, binomial(i+2,2)*B[i+1,j-1])  (second term always added)
    if i == 0:
        r = Fraction(int(j == 1))
    else:
        r = Fraction((i + j) * (6 * j + 2 * i - 5) * j * (2 * i + j - 1), 2 * (2 * i + 1) * (2 * j + i - 2)) * B(i, j)
    if j > 1:
        r += comb(i + 2, 2) * B(i + 1, j - 1)
    return r


def neglog1m(A):
    """-log(1 - A) for A with zero constant term."""
    N = A.N
    out = BS(N)
    p = BS.const(N, 1)
    for n in range(1, N + 1):
        p = p * A
        out = out + p * Fraction(1, n)
    return out


def rot_big(k, tau, N, nets=None):
    """c_tau(v', f') for k >= 3 as a BS in (X, Y) = quotient vertex/face orbit counts."""
    D, S, P, H = nets or networks(N)
    pv, pf = {'VF': (1, 1), 'VV': (2, 0), 'FF': (0, 2)}[tau]
    L = BS(N)
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
            e = k * (a + b)
            L[a, b] = phi(k) * Fraction(q) / (4 * e)
    u = BS.mono(N, 1, 0)
    w = BS.mono(N, 0, 1)
    if tau == 'FF':
        L = L - neglog1m(u * (1 + P + H)) * Fraction(phi(k), 4 * k)
    if tau == 'VV':
        L = L - neglog1m(w * (1 + S + H)) * Fraction(phi(k), 4 * k)
    return invert_xy(L, u * D, w * D), (pv, pf)


def main():
    import json
    import sys
    N = 12
    nets = networks(N)
    stats = {}
    for fn in sys.argv[1:]:
        for key, v in json.load(open(fn)).items():
            stats.setdefault(key, {}).update({t: Fraction(w) for t, w in v.items()})
    vmax = max(int(k.split(',')[0]) for k in stats)
    bad = 0
    checked = 0
    for k in range(3, 12):
        for tau in ('VF', 'VV', 'FF'):
            c, (pv, pf) = rot_big(k, tau, N, nets)
            pol = ''.join(sorted(tau))
            # compare on all (v,f) with v <= vmax (brute data complete there)
            for key, d in stats.items():
                v, f = map(int, key.split(','))
                want = d.get(json.dumps(['rot', k, pol]), 0)
                got = 0
                if (v - pv) % k == 0 and (f - pf) % k == 0:
                    a, b = (v - pv) // k, (f - pf) // k
                    got = c[a, b] if a + b <= N else None
                if got is None:
                    continue
                checked += 1
                if got != want:
                    bad += 1
                    print('MISMATCH', k, tau, v, f, got, want)
    print('checked', checked, 'bad', bad)


if __name__ == '__main__':
    main()
