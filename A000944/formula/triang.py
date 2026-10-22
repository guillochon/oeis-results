"""Symmetric triangulations of the sphere (simplicial polyhedra) by vertices, for A000109.

A000109(n) = (A253882(n) + achiral(n)) / 2, achiral = 2 (W_ref + W_rr) on the triangulation diagonal
f = 2v - 4 (U = S/2 + W_ref + W_rr). Here every series is in one variable t (vertices of M), so the
cost is far below the full vertex-by-face box of the A000944 run.

Reflections. In a simple triangulation M with a reflection, every point of the mirror circle is a vertex
(label V) or the midpoint of a crossing edge (a, sigma a); the end a on one side is labelled X. Every
face the mirror enters is mirror-symmetric, so the mirror runs V-V along an edge or V-X-V through the
two triangles on a crossing edge. The half H is a simple triangulated disk whose boundary is the cyclic
V/X word (no plain boundary vertices), with
    (1) no two X adjacent on the boundary,  (2) no interior chord joining two V
(a V-V chord and its mirror image would be a double edge; nothing else can make M non-simple, and a
simple triangulation with >= 4 vertices is 3-connected). M has #V + 2 #X + 2 #interior vertices.
H splits along its chords into chordless pieces: C(m, i) rooted chordless triangulations of an m-gon
with i interior vertices, by inverting Brown's count of rooted triangulations of a polygon. Arcs between
boundary labels a, b:  E = J + K o sum_m C_m E (W E)^(m-2), W = diag(alpha t, alpha t^2) (J: no XX
edge, K: no VV chord), and the root piece gives R = sum_m C_m tr(J (W E)^(m-1) W), rooted at a boundary
edge in the plane orientation, alpha marking the boundary length A. W_ref = sum_A [alpha^A] R / (4A).
"""
import sys

import numpy as np

from fastmul import mulmod
from t3 import PRIMES24


def fact_tables(M, p):
    f = [1] * (M + 1)
    for i in range(1, M + 1):
        f[i] = f[i - 1] * i % p
    fi = [pow(x, p - 2, p) for x in f]
    return f, fi


def brown(N, K, p):
    """psi[k][i]: rooted simple triangulations of a k-gon (chords allowed) with i interior vertices,
    k = 3..K, i = 0..N. Brown (1964): psi(n, m) = 2 (2m+3)! (4n+2m+1)! / ((m+2)! m! n! (3n+2m+3)!), k = m+3."""
    f, fi = fact_tables(4 * N + 2 * K + 8, p)
    psi = np.zeros((K + 1, N + 1), dtype=np.int64)
    for k in range(3, K + 1):
        m = k - 3
        for n in range(N + 1):
            psi[k, n] = (2 * f[2 * m + 3] * f[4 * n + 2 * m + 1] % p * fi[m + 2] % p * fi[m] % p * fi[n] % p
                         * fi[3 * n + 2 * m + 3]) % p
    return psi


def chordless(N, K, p):
    """C[m][i]: rooted chordless triangulations of an m-gon with i interior vertices (m = 3..K).
    T = sum_m C_m z A^(m-1), A = z + T/z, with z per boundary edge (rooted at a boundary edge)."""
    psi = brown(N, K, p)
    sh = (N + 1,)
    # B(z) = 1 + T/z^2 as rows by z-degree: B[d] = psi[d + 2] for d >= 1
    B = np.zeros((K + 1, N + 1), dtype=np.int64)
    B[0, 0] = 1
    for d in range(1, K - 1):
        B[d] = psi[d + 2]

    def zmul(X, Y):            # product of series in z (rows) with coefficients series in u
        out = np.zeros_like(X)
        for a in range(K + 1):
            if not X[a].any():
                continue
            for b in range(K + 1 - a):
                if Y[b].any():
                    out[a + b] = (out[a + b] + mulmod(X[a], Y[b], p, sh)) % p
        return out
    # [z^(k-1)] z^(m-1) B^(m-1) = [z^(k-m)] B^(m-1)
    Bpow = [np.zeros_like(B)]
    Bpow[0][0, 0] = 1
    for _ in range(K):
        Bpow.append(zmul(Bpow[-1], B))
    C = np.zeros((K + 1, N + 1), dtype=np.int64)
    for k in range(3, K + 1):
        acc = psi[k].copy()
        for m in range(3, k):
            acc = (acc - mulmod(C[m], Bpow[m - 1][k - m], p, sh)) % p
        C[k] = acc
    return C


def wref_triang(N, p):
    """W_ref(n) mod p for simple triangulations with n <= N vertices: array [n]."""
    K = N + 2                                   # boundary length <= number of M vertices
    I = N // 2 + 1                              # interior vertices of a piece: 2 i <= N
    C = chordless(I, K, p)
    sh = (N + 1, N + 1)                         # [t-degree, alpha-degree]

    def ser_t2(c):                              # sum_i c[i] t^(2i), as a [t, alpha] array
        a = np.zeros(sh, dtype=np.int64)
        for i, x in enumerate(c):
            if 2 * i <= N:
                a[2 * i, 0] = x
        return a
    Cm = {m: ser_t2(C[m]) for m in range(3, K + 1) if C[m].any()}
    mmax = max(Cm)
    zero = np.zeros(sh, dtype=np.int64)
    one = zero.copy(); one[0, 0] = 1
    wV = zero.copy(); wV[1, 1] = 1              # alpha t
    wX = zero.copy(); wX[2, 1] = 1              # alpha t^2
    Wd = [wV, wX]
    J = [[1, 1], [1, 0]]                        # boundary edge allowed (no X-X)
    Kc = [[0, 1], [1, 1]]                       # chord allowed (no V-V)

    def mm(A, B):
        return [[(mulmod(A[i][0], B[0][j], p, sh) + mulmod(A[i][1], B[1][j], p, sh)) % p for j in range(2)]
                for i in range(2)]

    def WE_of(E):
        return [[mulmod(Wd[i], E[i][j], p, sh) for j in range(2)] for i in range(2)]

    def horner(X, lo):
        """sum_{m >= lo} C_m X^(m - lo) as a 2x2 matrix of series."""
        acc = [[Cm.get(mmax, zero) if i == j else zero for j in range(2)] for i in range(2)]
        for m in range(mmax - 1, lo - 1, -1):
            acc = mm(acc, X)
            c = Cm.get(m, zero)
            acc = [[(acc[i][j] + (c if i == j else 0)) % p for j in range(2)] for i in range(2)]
        return acc
    E = [[one * J[i][j] for j in range(2)] for i in range(2)]
    for it in range(N + 3):
        WE = WE_of(E)
        S = mm(E, horner(WE, 3))                # sum_m C_m E (W E)^(m-3) ... times one more WE below
        S = mm(S, WE)                           # = sum_m C_m E (W E)^(m-2)
        new = [[(one * J[i][j] + S[i][j] * Kc[i][j]) % p for j in range(2)] for i in range(2)]
        if all(np.array_equal(new[i][j], E[i][j]) for i in range(2) for j in range(2)):
            break
        E = new
    else:
        raise RuntimeError('arc fixed point did not converge')
    WE = WE_of(E)
    T = mm(horner(WE, 3), WE)                   # sum_m C_m (W E)^(m-2) ... times WE twice more
    T = mm(T, WE)                               # sum_m C_m (W E)^(m-1)
    Wm = [[Wd[0], zero], [zero, Wd[1]]]
    T = mm(T, Wm)                               # (W E)^(m-1) W
    R = zero.copy()
    for i in range(2):
        for j in range(2):
            if J[j][i]:
                R = (R + T[i][j]) % p           # tr(J X) = sum_{i,j} J[j][i] X[i][j]
    out = np.zeros(N + 1, dtype=np.int64)
    for A in range(1, N + 1):
        col = R[:, A]
        if col.any():
            out = (out + col * pow(4 * A, p - 2, p)) % p   # rooted at an oriented boundary edge (one
            # orientation of the disk), so each half has A rootings; with the 1/2 of W_ref: 1/(4A)
    return out


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    p = PRIMES24[0]
    w = wref_triang(N, p)
    ref = np.load(f'wref_40_76_{p}.npy')
    ok = [n for n in range(4, min(N, 40) + 1) if int(w[n]) == int(ref[n, 2 * n - 4])]
    bad = [n for n in range(4, min(N, 40) + 1) if int(w[n]) != int(ref[n, 2 * n - 4])]
    print('ok', ok, '\nbad', bad)
    for n in bad[:6]:
        a, b = int(w[n]), int(ref[n, 2 * n - 4])
        print(n, 'ours', a, 'stored', b, 'ratio', a * pow(b, p - 2, p) % p if b else None)
