/* Strategy (b) of Mallows's 3D Apollonian packing (A154643, A154644, A154645).
 *
 * A sphere is an integer row vector v in Z^5 (sphere = v*C, C the starting
 * configuration in augmented curvature-center coordinates); its bend is
 * v0+v1+v2.  A quintuple is an ordered 5-tuple of spheres = the rows of a
 * group element; the row order is determined by the set, so it is a key.
 * Descartes reflection i: row i -> (sum of rows) - 2*row i.
 *
 * Level g is closed by reflecting every quintuple of generation g; a result
 * whose youngest sphere has generation <= g is another gen-g quintuple
 * ("extra"), one with a new sphere goes to level g+1.  Every quintuple ever
 * seen is kept, and a newly found quintuple of generation < g is counted in
 * `late` (it would mean the level-by-level closure is not exact).
 *
 * usage: stratb G   (closes levels 0..G-1: spheres to gen G, quintuples to gen G-1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef __int128 i128;
typedef struct { int32_t c[5]; } Sph;

/* ---- sphere table: open addressing on coordinates ---- */
static Sph *sph;          /* coordinates by id */
static uint8_t *sgen;     /* generation by id */
static uint32_t nsph, capsph;
static uint32_t *stab; static uint64_t smask;

static uint64_t hsph(const int64_t *c) {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (int k = 0; k < 5; k++) { h ^= (uint64_t)c[k] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); h *= 0xBF58476D1CE4E5B9ULL; }
    return h ^ (h >> 31);
}
static int sph_eq(const int32_t *a, const int64_t *c) { for (int k = 0; k < 5; k++) if (a[k] != c[k]) return 0; return 1; }
/* returns id; *isnew set */
static uint32_t sfind(const int64_t *c, int *isnew) {
    uint64_t i = hsph(c) & smask;
    for (;;) {
        uint32_t id = stab[i];
        if (id == UINT32_MAX) break;
        if (sph_eq(sph[id].c, c)) { *isnew = 0; return id; }
        i = (i + 1) & smask;
    }
    if (nsph >= capsph) { fprintf(stderr, "sphere capacity\n"); exit(1); }
    for (int k = 0; k < 5; k++) sph[nsph].c[k] = (int32_t)c[k];
    stab[i] = nsph; *isnew = 1;
    return nsph++;
}

/* ---- quintuple set: 5 sphere ids ---- */
typedef struct { uint32_t s[5]; } Q;
static Q *qtab; static uint64_t qmask, nq;
static uint64_t hq(const Q *q) {
    uint64_t h = 0xC2B2AE3D27D4EB4FULL;
    for (int k = 0; k < 5; k++) { h ^= q->s[k]; h *= 0x9E3779B97F4A7C15ULL; h ^= h >> 29; }
    return h;
}
static int qinsert(const Q *q) { /* 1 if new */
    uint64_t i = hq(q) & qmask;
    for (;;) {
        if (qtab[i].s[0] == UINT32_MAX) break;
        if (!memcmp(&qtab[i], q, 20)) return 0;
        i = (i + 1) & qmask;
    }
    if (++nq > (qmask + 1) / 10 * 7) { fprintf(stderr, "quint capacity\n"); exit(1); }
    qtab[i] = *q; return 1;
}

typedef struct { Q *a; uint64_t n, cap; } Vec;
static void push(Vec *v, const Q *q) {
    if (v->n == v->cap) { v->cap = v->cap ? v->cap * 2 : 1024; v->a = realloc(v->a, v->cap * sizeof(Q)); if (!v->a) { fprintf(stderr, "oom\n"); exit(1); } }
    v->a[v->n++] = *q;
}

static void pr128(i128 x) {
    char b[50]; int n = 0; int neg = x < 0; if (neg) x = -x;
    do { b[n++] = '0' + (int)(x % 10); x /= 10; } while (x);
    if (neg) putchar('-');
    while (n) putchar(b[--n]);
}

int main(int argc, char **argv) {
    int G = argc > 1 ? atoi(argv[1]) : 8;
    capsph = argc > 2 ? (uint32_t)strtoul(argv[2], 0, 10) : 1u << 24;
    uint64_t qcap = argc > 3 ? strtoull(argv[3], 0, 10) : 1ULL << 24;
    uint64_t sts = 1; while (sts < (uint64_t)capsph * 10 / 7) sts <<= 1;
    smask = sts - 1;
    stab = malloc(sts * 4); memset(stab, 0xff, sts * 4);
    sph = malloc((uint64_t)capsph * sizeof(Sph)); sgen = malloc(capsph);
    uint64_t qts = 1; while (qts < qcap * 10 / 7) qts <<= 1;
    qmask = qts - 1; qtab = malloc(qts * sizeof(Q)); memset(qtab, 0xff, qts * sizeof(Q));
    if (!stab || !sph || !sgen || !qtab) { fprintf(stderr, "alloc\n"); return 1; }

    Q q0; int isnew;
    for (int j = 0; j < 5; j++) { int64_t c[5] = {0}; c[j] = 1; q0.s[j] = sfind(c, &isnew); sgen[q0.s[j]] = 0; }
    qinsert(&q0);
    Vec cur = {0}, nxt = {0};
    push(&cur, &q0);
    uint64_t *cntS = calloc(G + 2, 8), *cntQ = calloc(G + 2, 8), late = 0;
    i128 *sumS = calloc(G + 2, 16), *sumT = calloc(G + 2, 16);
    cntS[0] = 5; sumS[0] = 3; sumT[0] = 3;
    int64_t maxabs = 0;

    for (int g = 0; g < G; g++) {
        nxt.n = 0;
        /* cur grows while we scan it: extra quintuples of generation g are appended */
        for (uint64_t qi = 0; qi < cur.n; qi++) {
            Q q = cur.a[qi];
            int64_t tot[5] = {0};
            for (int j = 0; j < 5; j++) for (int k = 0; k < 5; k++) tot[k] += sph[q.s[j]].c[k];
            for (int i = 0; i < 5; i++) {
                int64_t c[5];
                for (int k = 0; k < 5; k++) {
                    c[k] = tot[k] - 2 * sph[q.s[i]].c[k];
                    int64_t a = c[k] < 0 ? -c[k] : c[k]; if (a > maxabs) maxabs = a;
                }
                uint32_t id = sfind(c, &isnew);
                if (isnew) {
                    sgen[id] = g + 1;
                    int64_t b = c[0] + c[1] + c[2];
                    cntS[g + 1]++; sumS[g + 1] += b; sumT[g + 1] += (i128)b * b;
                }
                Q q2 = q; q2.s[i] = id;
                if (!qinsert(&q2)) continue;
                int gq = 0;
                for (int j = 0; j < 5; j++) if (sgen[q2.s[j]] > gq) gq = sgen[q2.s[j]];
                if (gq < g) late++;
                else if (gq == g) push(&cur, &q2);
                else push(&nxt, &q2);
            }
        }
        cntQ[g] = cur.n;
        Vec t = cur; cur = nxt; nxt = t;
        fprintf(stderr, "level %d closed: quints %llu, next %llu, spheres %u, maxabs %lld, late %llu\n",
                g, (unsigned long long)cntQ[g], (unsigned long long)cur.n, nsph, (long long)maxabs, (unsigned long long)late);
        if (maxabs > (1LL << 30)) { fprintf(stderr, "overflow risk\n"); return 1; }
    }
    printf("{\"G\": %d, \"late\": %llu,\n \"quintuples_b\": [", G, (unsigned long long)late);
    for (int g = 0; g < G; g++) printf("%s%llu", g ? ", " : "", (unsigned long long)cntQ[g]);
    printf("],\n \"A154643\": [");
    for (int g = 0; g <= G; g++) printf("%s%llu", g ? ", " : "", (unsigned long long)cntS[g]);
    printf("],\n \"A154644\": [");
    for (int g = 0; g <= G; g++) { if (g) printf(", "); if (sumS[g] % 3) { printf("\"x\""); continue; } pr128(sumS[g] / 3); }
    printf("],\n \"A154645\": [");
    for (int g = 0; g <= G; g++) { if (g) printf(", "); if (sumT[g] % 3) { printf("\"x\""); continue; } pr128(sumT[g] / 3); }
    printf("]}\n");
    return 0;
}
