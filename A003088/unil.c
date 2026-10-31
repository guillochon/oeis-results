/* A003088 (unilateral digraphs) and A035512 (strong digraphs) modulo one prime: the same method as
 * unilateral.py, in C, so that CRT over several primes gives the exact terms to larger n.
 *
 *   unil N P  ->  lines "n strong(n) mod P  a(n) mod P" for n = 0..N
 *
 * Partitions of 0..N are indexed by size; a partition is looked up by an additive 64-bit hash
 * (sum of a random word per part), so the hash of a multiset union is the sum of the hashes.
 * g(l,m) = sum gcd(a,b); h[l][b] = sum_{a in l} gcd(a,b), so g(l,m) = sum_{b in m} h[l][b].
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

static int N;
static u64 P;
static int np;                 /* number of partitions of 0..N */
static int *start;             /* start[s] = first index of size s; start[N+1] = np */
static int *psize, *plen;      /* size, number of parts */
static int *dpoff;             /* distinct parts of partition i: dpart/dmult[dpoff[i] .. dpoff[i+1]) */
static int *dpart, *dmult;
static u64 *H;                 /* additive hash */
static u64 R[1024];
static int *gg;                /* g(l,l) */
static unsigned short *hv;     /* hv[i*(N+1)+b] */
static u64 *htkey; static int *htval; static u64 htmask;

static u64 mulm(u64 a, u64 b) { return (u64)((u128)a * b % P); }
static u64 addm(u64 a, u64 b) { u64 c = a + b; return c >= P ? c - P : c; }
static u64 subm(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static u64 powm(u64 a, u64 e) { u64 r = 1; a %= P; while (e) { if (e & 1) r = mulm(r, a); a = mulm(a, a); e >>= 1; } return r; }
static u64 invm(u64 a) { return powm(a, P - 2); }

static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

static int look(u64 key) {
    u64 i = (key * 0x9E3779B97F4A7C15ULL) >> 20 & htmask;
    while (htval[i] >= 0) { if (htkey[i] == key) return htval[i]; i = (i + 1) & htmask; }
    fprintf(stderr, "hash miss\n"); exit(1);
}
static void ins(u64 key, int v) {
    u64 i = (key * 0x9E3779B97F4A7C15ULL) >> 20 & htmask;
    while (htval[i] >= 0) { if (htkey[i] == key) { fprintf(stderr, "hash collision\n"); exit(1); } i = (i + 1) & htmask; }
    htkey[i] = key; htval[i] = v;
}

/* partition enumeration (parts descending) */
static int cap, ndp, dcap;
static int cur[1024];
static void emit(int len, int size) {
    if (np == cap) {
        cap *= 2;
        psize = realloc(psize, cap * sizeof *psize); plen = realloc(plen, cap * sizeof *plen);
        dpoff = realloc(dpoff, (cap + 1) * sizeof *dpoff); H = realloc(H, cap * sizeof *H);
    }
    psize[np] = size; plen[np] = len; dpoff[np] = ndp;
    u64 h = 0;
    for (int i = 0; i < len; i++) h += R[cur[i]];
    H[np] = h;
    for (int i = 0; i < len;) {
        int j = i; while (j < len && cur[j] == cur[i]) j++;
        if (ndp + 1 >= dcap) { dcap *= 2; dpart = realloc(dpart, dcap * sizeof *dpart); dmult = realloc(dmult, dcap * sizeof *dmult); }
        dpart[ndp] = cur[i]; dmult[ndp] = j - i; ndp++;
        i = j;
    }
    np++;
}
static void gen(int rem, int maxp, int len, int size) {
    if (rem == 0) { emit(len, size); return; }
    for (int k = rem < maxp ? rem : maxp; k >= 1; k--) { cur[len] = k; gen(rem - k, k, len + 1, size); }
}

static int dot(const unsigned short *h, int mu) {
    int s = 0;
    for (int t = dpoff[mu]; t < dpoff[mu + 1]; t++) s += h[dpart[t]] * dmult[t];
    return s;
}

/* E step: enumerate sub-multisets of nu */
static u64 *E, *D, *pw;
static int snu, sdp0, sdp1;
static u64 sacc;
static int sub_l[64];
static void subrec(int t, u64 hl, u64 hm, int any_m) {
    if (t == sdp1) {
        if (!any_m) return;
        int l = look(hl), m = look(hm);
        int g = (gg[snu] - gg[l] - gg[m]) / 2;
        sacc = addm(sacc, mulm(mulm(E[l], D[m]), pw[g]));
        return;
    }
    int part = dpart[t], c = dmult[t];
    for (int take = 0; take <= c; take++)
        subrec(t + 1, hl + (u64)take * R[part], hm + (u64)(c - take) * R[part], any_m || take < c);
}

/* the top layers N-KH .. N-1: open-addressing maps (reduced Q index << 32 | reduced lam index) -> value */
#define KH 3
static u64 Lr[KH + 1][1024];
static int *redp[KH + 1];      /* redp[r][i]: index of partition i with its parts reduced by L_r */
static int *cid[KH + 1], *rid[KH + 1];  /* compact id of redp[r][i], and back */
static u64 *accr[KH + 1];
static int *touched, ntouch, ntouchmax;
static u64 *hl_q[1024], *hl_val[1024];  /* key: (additive hash of reduced Q, index of reduced lam >= 1) */
static int *hl_l[1024];
static size_t hl_cap[1024], hl_n[1024];
static u64 *rh[KH + 1];        /* rh[r][i]: hash of partition i with its parts reduced by L_r */

static u64 redhash(int i, int r) {
    u64 h = 0;
    for (int t = dpoff[i]; t < dpoff[i + 1]; t++) h += (u64)dmult[t] * R[Lr[r][dpart[t]]];
    return h;
}

static void hadd(int s, u64 q, int l, u64 w);
static void hgrow(int s) {
    size_t oc = hl_cap[s], nc = oc ? 2 * oc : 1024;
    u64 *oq = hl_q[s], *ov = hl_val[s];
    int *ol = hl_l[s];
    hl_q[s] = malloc(nc * sizeof(u64)); hl_val[s] = calloc(nc, sizeof(u64)); hl_l[s] = calloc(nc, sizeof(int));
    hl_cap[s] = nc; hl_n[s] = 0;
    for (size_t i = 0; i < oc; i++) if (ol[i]) hadd(s, oq[i], ol[i], ov[i]);
    free(oq); free(ov); free(ol);
}
static void hadd(int s, u64 q, int l, u64 w) {
    if (2 * (hl_n[s] + 1) > hl_cap[s]) hgrow(s);
    size_t mask = hl_cap[s] - 1, i = (size_t)(((q + (u64)l * 0x9E3779B97F4A7C15ULL) * 0xD6E8FEB86659FD93ULL) >> 17) & mask;
    while (hl_l[s][i] && (hl_l[s][i] != l || hl_q[s][i] != q)) i = (i + 1) & mask;
    if (!hl_l[s][i]) { hl_q[s][i] = q; hl_l[s][i] = l; hl_n[s]++; }
    hl_val[s][i] = addm(hl_val[s][i], w);
}

static int mobius(int n) {
    int r = 1;
    for (int p = 2; p * p <= n; p++) if (n % p == 0) { n /= p; if (n % p == 0) return 0; r = -r; }
    return n > 1 ? -r : r;
}

int main(int argc, char **argv) {
    N = atoi(argv[1]);
    P = strtoull(argv[2], 0, 10);
    clock_t t0 = clock();
    u64 seed = 0x243F6A8885A308D3ULL;
    for (int i = 0; i < 1024; i++) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; R[i] = seed; }

    cap = 1024; dcap = 4096;
    psize = malloc(cap * sizeof *psize); plen = malloc(cap * sizeof *plen);
    dpoff = malloc((cap + 1) * sizeof *dpoff); H = malloc(cap * sizeof *H);
    dpart = malloc(dcap * sizeof *dpart); dmult = malloc(dcap * sizeof *dmult);
    start = malloc((N + 2) * sizeof *start);
    for (int s = 0; s <= N; s++) { start[s] = np; gen(s, s, 0, s); }
    start[N + 1] = np; dpoff[np] = ndp;

    u64 hs = 1; while (hs < 4 * (u64)np) hs <<= 1;
    htmask = hs - 1; htkey = malloc(hs * sizeof *htkey); htval = malloc(hs * sizeof *htval);
    for (u64 i = 0; i < hs; i++) htval[i] = -1;
    for (int i = 0; i < np; i++) ins(H[i], i);

    hv = calloc((size_t)np * (N + 1), sizeof *hv);
    gg = malloc(np * sizeof *gg);
    for (int i = 0; i < np; i++) {
        unsigned short *h = hv + (size_t)i * (N + 1);
        for (int b = 1; b <= N; b++) {
            int s = 0;
            for (int t = dpoff[i]; t < dpoff[i + 1]; t++) s += gcd(dpart[t], b) * dmult[t];
            h[b] = s;
        }
        gg[i] = dot(h, i);
    }

    pw = malloc((N * N + 2) * sizeof *pw);
    pw[0] = 1; for (int e = 1; e <= N * N + 1; e++) pw[e] = addm(pw[e - 1], pw[e - 1]);

    /* 1. digraphs: D = 2^(g(l,l) - len) / z */
    D = malloc(np * sizeof *D);
    for (int i = 0; i < np; i++) {
        u64 z = 1;
        for (int t = dpoff[i]; t < dpoff[i + 1]; t++) {
            int k = dpart[t], c = dmult[t];
            for (int j = 1; j <= c; j++) z = mulm(z, (u64)k * j);
        }
        D[i] = mulm(pw[gg[i] - plen[i]], invm(z));
    }
    /* 2. E (.) D = 1 */
    E = malloc(np * sizeof *E);
    E[0] = 1;
    for (int i = 1; i < np; i++) {
        snu = i; sdp0 = dpoff[i]; sdp1 = dpoff[i + 1]; sacc = 0;
        subrec(sdp0, 0, 0, 0);
        E[i] = subm(0, sacc);
    }
    /* 3. Lg = log E, L = -Lg */
    u64 *Lg = calloc(np, sizeof *Lg), *acc = calloc(np, sizeof *acc);
    for (int n = 1; n <= N; n++) {
        for (int i = start[n]; i < start[n + 1]; i++) acc[i] = mulm(n, E[i]);
        for (int k = 1; k < n; k++)
            for (int a = start[k]; a < start[k + 1]; a++) {
                if (!Lg[a]) continue;
                u64 la = mulm(k, Lg[a]);
                for (int b = start[n - k]; b < start[n - k + 1]; b++) {
                    int nu = look(H[a] + H[b]);
                    acc[nu] = subm(acc[nu], mulm(la, E[b]));
                }
            }
        u64 ni = invm(n);
        for (int i = start[n]; i < start[n + 1]; i++) Lg[i] = mulm(acc[i], ni);
    }
    /* Z_S = sum_m mu(m)/m L[p_i -> p_im] */
    u64 *Z = calloc(np, sizeof *Z);
    for (int m = 1; m <= N; m++) {
        int mu = mobius(m);
        if (!mu) continue;
        u64 w = invm(m); if (mu < 0) w = subm(0, w);
        w = subm(0, w);  /* L = -Lg */
        for (int i = 1; i < start[N / m + 1]; i++) {
            u64 h = 0;
            for (int t = dpoff[i]; t < dpoff[i + 1]; t++) h += (u64)dmult[t] * R[m * dpart[t]];
            int j = look(h);
            Z[j] = addm(Z[j], mulm(w, Lg[i]));
        }
    }
    free(Lg); free(acc); free(E); free(D);
    u64 *strong = calloc(N + 1, sizeof *strong);
    strong[0] = 1;
    for (int n = 1; n <= N; n++) for (int i = start[n]; i < start[n + 1]; i++) strong[n] = addm(strong[n], Z[i]);
    fprintf(stderr, "# N=%d cycle index %.1fs, %d partitions\n", N, (double)(clock() - t0) / CLOCKS_PER_SEC, np);

    /* 4. chain DP over states (Q, lam) in layers s = |Q| + |lam|.
     * Layers s <= M = N - 1 - KH: dense, F[off[q] + lam - 1] for lam in 1..start[M-|Q|+1]-1.
     * Layers N-KH .. N-1: what is left to append has size r = N - s <= KH, and a state only matters through
     *   gcd(a, b) for its parts a and b <= r, i.e. through the parts reduced to L_r(a) = lcm{d <= r : d | a};
     *   these layers are hash maps keyed by (reduced Q, reduced lam).
     * Layer N: never stored, only summed. */
    int M = N - 1 - KH;
    size_t *off = malloc((np + 1) * sizeof *off);
    size_t tot = 0;
    for (int q = 0; q < np; q++) { off[q] = tot; if (psize[q] < M) tot += (size_t)start[M - psize[q] + 1] - 1; }
    off[np] = tot;
    u64 *F = calloc(tot ? tot : 1, sizeof *F);
    if (!F) { fprintf(stderr, "out of memory (%zu states)\n", tot); return 1; }
    for (int r = 1; r <= KH; r++) {  /* L_r and the reduced index of every partition */
        for (int a = 1; a <= N; a++) {
            u64 L = 1;
            for (int d = 1; d <= r && d <= a; d++) if (a % d == 0) L = L / gcd(L, d) * d;
            Lr[r][a] = L;
        }
        redp[r] = malloc(np * sizeof(int)); rh[r] = malloc(np * sizeof(u64));
        for (int i = 0; i < np; i++) { rh[r][i] = redhash(i, r); redp[r][i] = look(rh[r][i]); }
        /* compact ids of the reduced partitions, for the per-state accumulators */
        int *map = malloc(np * sizeof(int)), k = 0;
        for (int i = 0; i < np; i++) map[i] = -1;
        cid[r] = malloc(np * sizeof(int)); rid[r] = malloc(np * sizeof(int));
        for (int i = 0; i < np; i++) {
            int j = redp[r][i];
            if (map[j] < 0) { map[j] = k; rid[r][k++] = j; }
            cid[r][i] = map[j];
        }
        accr[r] = calloc(k, sizeof(u64));
        ntouchmax += k;
        free(map);
    }
    touched = malloc((ntouchmax + 1) * sizeof(int));
    u64 *total = calloc(N + 1, sizeof *total);
    total[0] = 1;
    for (int l = 1; l < np; l++) {
        int t = psize[l];
        if (t <= M) F[off[0] + l - 1] = Z[l];
        else if (t < N) hadd(t, 0, redp[N - t][l], Z[l]);  /* reduced empty Q has hash 0 */
        else total[N] = addm(total[N], Z[l]);
    }
    unsigned short *h2 = malloc((N + 1) * sizeof *h2);
    size_t nhash = 0;
    for (int s = 1; s < N; s++) {
        int r = N - s;
        int mlim = start[r + 1];          /* mu indices 1 .. mlim-1 */
        int mdense = start[M - s + 1 > 0 ? M - s + 1 : 0];  /* mu with s + |mu| <= M: indices < mdense */
        if (mdense < 1) mdense = 1;
        int mtop = start[r];              /* mu with s + |mu| = N: indices >= mtop */
        int dense = s <= M;
        /* the states of layer s: dense (q, l) pairs, or the entries of the hash map */
        size_t nst = dense ? (size_t)start[s] : hl_cap[s];
        for (size_t qi = 0; qi < nst; qi++) {
            int q, lbeg, lend;
            if (dense) { q = (int)qi; int ls = s - psize[q]; lbeg = start[ls]; lend = start[ls + 1]; }
            else { if (!hl_l[s][qi]) continue; q = look(hl_q[s][qi]); lbeg = hl_l[s][qi]; lend = lbeg + 1; }
            const unsigned short *hq = hv + (size_t)q * (N + 1);
            for (int l = lbeg; l < lend; l++) {
                u64 v = dense ? F[off[q] + l - 1] : hl_val[s][qi];
                if (!v) continue;
                total[s] = addm(total[s], v);
                const unsigned short *hl = hv + (size_t)l * (N + 1);
                for (int b = 1; b <= r; b++) h2[b] = hq[b] + hl[b];
                u64 hql = H[q] + H[l];
                int pn = (dense && s < M) ? look(hql) : -1;
                u64 pr[KH + 1];               /* hash of reduced Q + lam, for the hash layers */
                for (int rr = 1; rr <= KH && rr < r; rr++) pr[rr] = rh[rr][q] + rh[rr][l];
                u64 *Fp = pn >= 0 ? F + off[pn] - 1 : 0;
                for (int m = 1; m < mdense; m++) {
                    int gq = 0, g2 = 0;
                    for (int t = dpoff[m]; t < dpoff[m + 1]; t++) { gq += hq[dpart[t]] * dmult[t]; g2 += h2[dpart[t]] * dmult[t]; }
                    Fp[m] = addm(Fp[m], mulm(mulm(v, Z[m]), subm(pw[g2], pw[gq])));
                }
                /* hash layers: many mu share a reduced key, so sum Z[mu] * (2^g2 - 2^gq) per key, then one write */
                for (int m = mdense; m < mtop; m++) {
                    int gq = 0, g2 = 0;
                    for (int t = dpoff[m]; t < dpoff[m + 1]; t++) { gq += hq[dpart[t]] * dmult[t]; g2 += h2[dpart[t]] * dmult[t]; }
                    int rr = r - psize[m], c = cid[rr][m];
                    if (!accr[rr][c]) touched[ntouch++] = rr << 24 | c;
                    accr[rr][c] = addm(accr[rr][c], mulm(Z[m], subm(pw[g2], pw[gq])));
                    if (!accr[rr][c]) accr[rr][c] = P;  /* keep "touched" distinct from 0; P == 0 mod P */
                }
                for (int k = 0; k < ntouch; k++) {
                    int rr = touched[k] >> 24, c = touched[k] & 0xffffff;
                    u64 a = accr[rr][c] == P ? 0 : accr[rr][c];
                    accr[rr][c] = 0;
                    if (a) hadd(N - rr, pr[rr], rid[rr][c], mulm(v, a));
                }
                ntouch = 0;
                u64 acc = 0;                  /* layer N: only the sum */
                for (int m = mtop; m < mlim; m++) {
                    int gq = 0, g2 = 0;
                    for (int t = dpoff[m]; t < dpoff[m + 1]; t++) { gq += hq[dpart[t]] * dmult[t]; g2 += h2[dpart[t]] * dmult[t]; }
                    acc = addm(acc, mulm(Z[m], subm(pw[g2], pw[gq])));
                }
                total[N] = addm(total[N], mulm(v, acc));
            }
        }
        if (!dense) { nhash += hl_n[s]; free(hl_q[s]); free(hl_val[s]); free(hl_l[s]); hl_l[s] = 0; }
        if (s == M) { free(F); F = 0; }
        fprintf(stderr, "# s=%d %.1fs\n", s, (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    for (int n = 0; n <= N; n++) printf("%d %llu %llu\n", n, (unsigned long long)strong[n], (unsigned long long)total[n]);
    fprintf(stderr, "# N=%d done %.1fs, %zu dense states, %zu hashed states\n", N, (double)(clock() - t0) / CLOCKS_PER_SEC, tot, nhash);
    return 0;
}
