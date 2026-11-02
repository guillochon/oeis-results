/* Strategy (b), modulo the S5 symmetry permuting the five generators.
 *
 * Conjugating by a permutation matrix P (R_pi(i) = P R_i P^-1) maps group
 * elements to group elements, rows to rows (sphere v -> v P^-1, i.e. its
 * coordinates permuted) and preserves strategy-(b) generations, which do not
 * depend on the bends.  So the BFS of stratb.c runs on orbits:
 *   sphere class   = sorted coordinate vector, orbit size 120 / prod(mult!)
 *   quintuple class = lexicographically least P M P^-1 over the 120 P,
 *                     orbit size 120 / #minimizers; stored as a 128-bit fingerprint.
 * Bends of Mallows's start are v0+v1+v2; summed over all arrangements of a
 * class's coordinates (uniformly, #arr of them):
 *   sum b   = #arr * 3 S1 / 5,     sum b^2 = #arr * 3 (S2 + S1^2) / 10,
 * with S1 = sum v_k, S2 = sum v_k^2.
 *
 * usage: stratb_sym G [sphere_cap] [quint_cap]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef __int128 i128;
typedef struct { int32_t m[25]; } QM;   /* quintuple matrix, row-major */

static int perm[120][5], nperm;
static void genperms(int *a, int k) {
    if (k == 5) { memcpy(perm[nperm++], a, sizeof(int) * 5); return; }
    for (int i = k; i < 5; i++) {
        int t = a[k]; a[k] = a[i]; a[i] = t;
        genperms(a, k + 1);
        t = a[k]; a[k] = a[i]; a[i] = t;
    }
}

/* ---- sphere classes ---- */
typedef struct { int32_t c[5]; } Sph;
static Sph *sph; static uint8_t *sgen; static uint32_t nsph, capsph;
static uint32_t *stab; static uint64_t smask;
static uint64_t mix(uint64_t h) { h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33; return h; }
static uint64_t hs(const int32_t *c) { uint64_t h = 0x1234567; for (int k = 0; k < 5; k++) h = mix(h ^ (uint32_t)c[k] ^ ((uint64_t)k << 40)); return h; }
static void sort5(int32_t *c) { for (int i = 1; i < 5; i++) { int32_t x = c[i]; int j = i - 1; while (j >= 0 && c[j] > x) { c[j + 1] = c[j]; j--; } c[j + 1] = x; } }
static uint32_t sfind(const int32_t *v, int *isnew) {
    int32_t c[5]; memcpy(c, v, 20); sort5(c);
    uint64_t i = hs(c) & smask;
    for (;;) {
        uint32_t id = stab[i];
        if (id == UINT32_MAX) break;
        if (!memcmp(sph[id].c, c, 20)) { *isnew = 0; return id; }
        i = (i + 1) & smask;
    }
    if (nsph >= capsph) { fprintf(stderr, "sphere capacity\n"); exit(1); }
    memcpy(sph[nsph].c, c, 20); stab[i] = nsph; *isnew = 1;
    return nsph++;
}

/* ---- quintuple fingerprints ---- */
typedef struct { uint64_t a, b; } FP;
static FP *qtab; static uint64_t qmask, nq;
static int qinsert(FP f) {
    if (f.a == 0 && f.b == 0) f.a = 1;
    uint64_t i = f.a & qmask;
    for (;;) {
        if (qtab[i].a == 0 && qtab[i].b == 0) break;
        if (qtab[i].a == f.a && qtab[i].b == f.b) return 0;
        i = (i + 1) & qmask;
    }
    if (++nq > (qmask + 1) / 10 * 7) { fprintf(stderr, "quint capacity\n"); exit(1); }
    qtab[i] = f; return 1;
}

/* canonical form: least P M P^-1, N[p[i]][p[j]] = M[i][j]; returns #minimizers */
static int canon(const QM *M, QM *out) {
    int best = -1, nmin = 0;
    QM N;
    for (int t = 0; t < 120; t++) {
        const int *p = perm[t];
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) N.m[p[i] * 5 + p[j]] = M->m[i * 5 + j];
        int c = -1;
        if (best >= 0) {
            c = 0;
            for (int k = 0; k < 25; k++) if (N.m[k] != out->m[k]) { c = N.m[k] < out->m[k] ? -1 : 1; break; }
        }
        if (c < 0) { *out = N; best = t; nmin = 1; }
        else if (c == 0) nmin++;
    }
    return nmin;
}
static FP fp(const QM *q) {
    uint64_t a = 0x9E3779B97F4A7C15ULL, b = 0x2545F4914F6CDD1DULL;
    for (int k = 0; k < 25; k++) { a = mix(a ^ (uint32_t)q->m[k] ^ ((uint64_t)k << 37)); b = mix(b + (uint32_t)q->m[k] * 0x100000001B3ULL + k); }
    FP f = {a, b}; return f;
}

typedef struct { QM *a; uint64_t n, cap; } Vec;
static void push(Vec *v, const QM *q) {
    if (v->n == v->cap) { v->cap = v->cap ? v->cap + v->cap / 2 : 1024; v->a = realloc(v->a, v->cap * sizeof(QM)); if (!v->a) { fprintf(stderr, "oom\n"); exit(1); } }
    v->a[v->n++] = *q;
}
static void pr128(i128 x) {
    char b[50]; int n = 0; int neg = x < 0; if (neg) x = -x;
    do { b[n++] = '0' + (int)(x % 10); x /= 10; } while (x);
    if (neg) putchar('-');
    while (n) putchar(b[--n]);
}
static int arrangements(const int32_t *c) { /* c sorted */
    int r = 120, run = 1;
    static const int f[6] = {1, 1, 2, 6, 24, 120};
    for (int k = 1; k <= 5; k++) {
        if (k < 5 && c[k] == c[k - 1]) run++;
        else { r /= f[run]; run = 1; }
    }
    return r;
}

int main(int argc, char **argv) {
    int G = argc > 1 ? atoi(argv[1]) : 8;
    capsph = argc > 2 ? (uint32_t)strtoul(argv[2], 0, 10) : 1u << 22;
    uint64_t qcap = argc > 3 ? strtoull(argv[3], 0, 10) : 1ULL << 22;
    int a0[5] = {0, 1, 2, 3, 4}; genperms(a0, 0);
    uint64_t t = 1; while (t < (uint64_t)capsph * 10 / 7) t <<= 1;
    smask = t - 1; stab = malloc(t * 4); memset(stab, 0xff, t * 4);
    sph = malloc((uint64_t)capsph * sizeof(Sph)); sgen = malloc(capsph);
    t = 1; while (t < qcap * 10 / 7) t <<= 1;
    qmask = t - 1; qtab = calloc(t, sizeof(FP));
    if (!stab || !sph || !sgen || !qtab) { fprintf(stderr, "alloc\n"); return 1; }

    uint64_t *cntS = calloc(G + 2, 8), *cntQ = calloc(G + 2, 8), *orbQ = calloc(G + 2, 8), late = 0;
    i128 *sumS = calloc(G + 2, 16), *sumT = calloc(G + 2, 16);
    QM q0; memset(&q0, 0, sizeof q0);
    int isnew;
    for (int j = 0; j < 5; j++) { q0.m[j * 5 + j] = 1; uint32_t id = sfind(&q0.m[j * 5], &isnew); sgen[id] = 0; }
    cntS[0] = 5; sumS[0] = 5; sumT[0] = 10;
    qinsert(fp(&q0));
    Vec cur = {0}, nxt = {0};
    push(&cur, &q0);
    int32_t maxabs = 0;

    for (int g = 0; g < G; g++) {
        nxt.n = 0;
        int last = (g == G - 1);
        for (uint64_t qi = 0; qi < cur.n; qi++) {
            QM q = cur.a[qi];
            QM c;
            int st = canon(&q, &c);
            cntQ[g] += 120 / st; orbQ[g]++;
            int64_t tot[5] = {0};
            for (int j = 0; j < 5; j++) for (int k = 0; k < 5; k++) tot[k] += q.m[j * 5 + k];
            uint8_t rg[5];
            for (int j = 0; j < 5; j++) { uint32_t id = sfind(&q.m[j * 5], &isnew); if (isnew) { fprintf(stderr, "missing sphere\n"); return 1; } rg[j] = sgen[id]; }
            for (int i = 0; i < 5; i++) {
                int32_t v[5];
                for (int k = 0; k < 5; k++) {
                    int64_t x = tot[k] - 2 * (int64_t)q.m[i * 5 + k];
                    if (x > 2000000000 || x < -2000000000) { fprintf(stderr, "overflow\n"); return 1; }
                    v[k] = (int32_t)x;
                    int32_t ax = v[k] < 0 ? -v[k] : v[k]; if (ax > maxabs) maxabs = ax;
                }
                uint32_t id = sfind(v, &isnew);
                if (isnew) {
                    sgen[id] = g + 1;
                    const int32_t *s = sph[id].c;
                    int arr = arrangements(s);
                    i128 S1 = 0, S2 = 0;
                    for (int k = 0; k < 5; k++) { S1 += s[k]; S2 += (i128)s[k] * s[k]; }
                    cntS[g + 1] += arr;
                    sumS[g + 1] += (i128)arr * S1;              /* bend sum / s0 = this / 5 */
                    sumT[g + 1] += (i128)arr * (S2 + S1 * S1);  /* square sum / t0 = this / 10 */
                }
                int gq = sgen[id];
                for (int j = 0; j < 5; j++) if (j != i && rg[j] > gq) gq = rg[j];
                if (gq < g) { /* old quintuple: must already be known */
                    QM q2 = q; memcpy(&q2.m[i * 5], v, 20);
                    QM c2; canon(&q2, &c2);
                    if (qinsert(fp(&c2))) late++;
                    continue;
                }
                if (gq > g && last) continue;   /* level G not needed */
                QM q2 = q; memcpy(&q2.m[i * 5], v, 20);
                QM c2; canon(&q2, &c2);
                if (!qinsert(fp(&c2))) continue;
                if (gq == g) push(&cur, &c2); else push(&nxt, &c2);
            }
        }
        Vec tv = cur; cur = nxt; nxt = tv;
        fprintf(stderr, "level %d closed: quints %llu (%llu orbits), next orbits %llu, sphere classes %u, maxabs %d, late %llu\n",
                g, (unsigned long long)cntQ[g], (unsigned long long)orbQ[g], (unsigned long long)cur.n, nsph, maxabs, (unsigned long long)late);
        fflush(stderr);
    }
    printf("{\"G\": %d, \"late\": %llu,\n \"quintuples_b\": [", G, (unsigned long long)late);
    for (int g = 0; g < G; g++) printf("%s%llu", g ? ", " : "", (unsigned long long)cntQ[g]);
    printf("],\n \"A154643\": [");
    for (int g = 0; g <= G; g++) printf("%s%llu", g ? ", " : "", (unsigned long long)cntS[g]);
    printf("],\n \"A154644\": [");
    for (int g = 0; g <= G; g++) { if (g) printf(", "); if (sumS[g] % 5) { printf("\"x\""); continue; } pr128(sumS[g] / 5); }
    printf("],\n \"A154645\": [");
    for (int g = 0; g <= G; g++) { if (g) printf(", "); if (sumT[g] % 10) { printf("\"x\""); continue; } pr128(sumT[g] / 10); }
    printf("]}\n");
    return 0;
}
