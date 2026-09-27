/*
 * skewcol.c -- INDEPENDENT CHECK of U(m) and U'(m) from skew.c / sequences.py.
 *
 * Instead of Polya counting over automorphism groups of connected structures
 * and an Euler transform, generate the coloured skew posets directly and
 * count canonical forms (nauty with vertex colours).
 *
 *   mode u : points coloured by block size s>=1            -> U(m)  (A109459)
 *   mode s : literal 2i coloured (p,q), literal 2i+1 (q,p)  -> U'(m) (A109458)
 *
 * Weight of an object = sum of block sizes.  Objects of weight w are made
 * from objects of weight w-s by adding a pair {y,~y} with y maximal and
 * {l : ~y < l} = B a consistent up-set of the parent, coloured with size s.
 *
 * Build: gcc -O2 -fopenmp -I$NAUTY skewcol.c $NAUTY/nauty.a -lm -o skewcol
 * Run:   ./skewcol u 9      ./skewcol s 8
 * Output: "U w count" lines; count(0)=1 is implicit.
 */
#define MAXN 64
#include "nauty.h"
#include "nautinv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define WMAX 9
#define NVMAX 18
#define NSHARD 1024

typedef struct { uint64_t w[7]; } Key;   /* 324 bits rows + 18*6 bits colours */

static void putbits(Key *k, int pos, int nbits, uint64_t v) {
    int wi = pos >> 6, off = pos & 63;
    k->w[wi] |= v << off;
    if (off + nbits > 64) k->w[wi + 1] |= v >> (64 - off);
}
static uint64_t getbits(const Key *k, int pos, int nbits) {
    int wi = pos >> 6, off = pos & 63;
    uint64_t v = k->w[wi] >> off;
    if (off + nbits > 64) v |= k->w[wi + 1] << (64 - off);
    return v & ((1ull << nbits) - 1);
}
static void pack(const uint32_t *rows, const int *col, int nv, Key *k) {
    memset(k, 0, sizeof(Key));
    for (int v = 0; v < nv; v++) putbits(k, 18 * v, 18, rows[v] & 0x3FFFFu);
    for (int v = 0; v < nv; v++) putbits(k, 324 + 6 * v, 6, (uint64_t)col[v]);
}
static void unpack(const Key *k, int nv, uint32_t *rows, int *col) {
    for (int v = 0; v < nv; v++) rows[v] = (uint32_t)getbits(k, 18 * v, 18);
    for (int v = 0; v < nv; v++) col[v] = (int)getbits(k, 324 + 6 * v, 6);
}
static uint64_t hash_key(const Key *k) {
    uint64_t h = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < 7; i++) { h ^= k->w[i]; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 31; }
    return h;
}

typedef struct { Key *list; uint32_t n, cap; uint32_t *tab; uint32_t tcap; omp_lock_t lock; } Shard;
static void shard_init(Shard *s) { s->n = 0; s->cap = 256; s->list = malloc(s->cap * sizeof(Key)); s->tcap = 1024; s->tab = calloc(s->tcap, 4); omp_init_lock(&s->lock); }
static void shard_grow(Shard *s) {
    uint32_t nc = s->tcap * 2; uint32_t *nt = calloc(nc, 4);
    for (uint32_t i = 0; i < s->tcap; i++) if (s->tab[i]) {
        uint32_t j = (uint32_t)(hash_key(&s->list[s->tab[i] - 1]) >> 20) & (nc - 1);
        while (nt[j]) j = (j + 1) & (nc - 1);
        nt[j] = s->tab[i];
    }
    free(s->tab); s->tab = nt; s->tcap = nc;
}
static int shard_insert(Shard *s, const Key *k, uint64_t h) {
    uint32_t j = (uint32_t)(h >> 20) & (s->tcap - 1);
    while (s->tab[j]) { if (!memcmp(&s->list[s->tab[j] - 1], k, sizeof(Key))) return 0; j = (j + 1) & (s->tcap - 1); }
    if (s->n == s->cap) { s->cap *= 2; s->list = realloc(s->list, s->cap * sizeof(Key)); }
    s->list[s->n] = *k; s->tab[j] = ++s->n;
    if (s->n * 2 > s->tcap) shard_grow(s);
    return 1;
}

/* colour codes: mode u -> s (1..9); mode s -> 10*p+q style index p*10+q (<100 needs 7 bits) -> use p*(WMAX+1)+q < 100? use 6 bits: p,q<=9 -> p*10+q up to 99 > 63. Use (p+q)*(p+q+1)/2 + q <= 54 for p+q<=9. */
static int colcode(int p, int q) { int s = p + q; return s * (s + 1) / 2 + q; }

/* canonical form with vertex colours */
static void canon(const uint32_t *rows, const int *col, int nv, uint32_t *crows, int *ccol) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orb[MAXN];
    DEFAULTOPTIONS_DIGRAPH(opt); statsblk st;
    opt.getcanon = TRUE; opt.defaultptn = FALSE;
    EMPTYGRAPH(g, 1, nv);
    for (int v = 0; v < nv; v++) for (int u = 0; u < nv; u++) if (rows[v] >> u & 1) ADDONEARC(g, v, u, 1);
    /* lab sorted by colour, ptn marks cell ends */
    int idx = 0;
    for (int c = 0; c < 64; c++) {
        int start = idx;
        for (int v = 0; v < nv; v++) if (col[v] == c) { lab[idx] = v; ptn[idx] = 1; idx++; }
        if (idx > start) ptn[idx - 1] = 0;
    }
    densenauty(g, lab, ptn, orb, &opt, &st, 1, nv, cg);
    for (int v = 0; v < nv; v++) {
        uint32_t m = 0; set *r = GRAPHROW(cg, v, 1);
        for (int u = 0; u < nv; u++) if (ISELEMENT(r, u)) m |= 1u << u;
        crows[v] = m; ccol[v] = col[lab[v]];
    }
}

static void upsets_rec(const uint32_t *up, const int *order, int nv, int pos, uint32_t B, uint32_t *out, uint32_t *cnt) {
    if (pos == nv) { out[(*cnt)++] = B; return; }
    int l = order[pos];
    upsets_rec(up, order, nv, pos + 1, B, out, cnt);
    if ((up[l] & ~B) == 0 && !(B >> (l ^ 1) & 1)) upsets_rec(up, order, nv, pos + 1, B | (1u << l), out, cnt);
}

/* canonical form -> up[] and colours with pairs relabelled to (2i,2i+1) */
static void form_to_up(const Key *k, int nv, uint32_t *up, int *col) {
    uint32_t rows[NVMAX]; int c0[NVMAX], newlab[NVMAX], next = 0;
    unpack(k, nv, rows, c0);
    for (int v = 0; v < nv; v++) newlab[v] = -1;
    for (int v = 0; v < nv; v++) if (newlab[v] < 0) {
        int partner = -1;
        for (int u = 0; u < nv; u++) if ((rows[v] >> u & 1) && (rows[u] >> v & 1)) { partner = u; break; }
        if (partner < 0) { fprintf(stderr, "no partner\n"); exit(1); }
        newlab[v] = 2 * next; newlab[partner] = 2 * next + 1; next++;
    }
    for (int v = 0; v < nv; v++) {
        uint32_t m = 0;
        for (int u = 0; u < nv; u++) if (rows[v] >> u & 1) m |= 1u << newlab[u];
        up[newlab[v]] = m & ~(1u << (newlab[v] ^ 1));
        col[newlab[v]] = c0[v];
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: skewcol u|s WMAX\n"); return 1; }
    int signed_mode = argv[1][0] == 's';
    int wmax = atoi(argv[2]); if (wmax < 1 || wmax > WMAX) return 1;
    Shard *lev[WMAX + 1];
    for (int w = 1; w <= wmax; w++) { lev[w] = malloc(NSHARD * sizeof(Shard)); for (int s = 0; s < NSHARD; s++) shard_init(&lev[w][s]); }
    /* new-pair colour options for size s: unsigned -> {s}; signed -> (p,q), p+q=s */
    for (int w = 1; w <= wmax; w++) {
        double t0 = omp_get_wtime();
        uint64_t total = 0;
        for (int pw = 0; pw < w; pw++) {          /* parent weight */
            int s = w - pw;                        /* new block size */
            int ncol = signed_mode ? s + 1 : 1; int colopt[WMAX + 2];
            for (int q = 0; q < ncol; q++) colopt[q] = signed_mode ? colcode(s - q, q) : s;
            if (pw == 0) {   /* empty parent: single point */
                for (int q = 0; q < ncol; q++) {
                    uint32_t rows[2] = {2u, 1u}, crows[NVMAX]; int col[2], ccol[NVMAX]; Key k;
                    col[0] = colopt[q]; col[1] = signed_mode ? colcode(q, s - q) : s;
                    canon(rows, col, 2, crows, ccol); pack(crows, ccol, 2, &k);
                    uint64_t h = hash_key(&k); shard_insert(&lev[w][h & (NSHARD - 1)], &k, h);
                }
                continue;
            }
            /* parents of weight pw may have any number of points r: infer from rows: nv = 2r; we store nv in the key? derive: count vertices with nonzero row (pair arcs make every row nonzero) */
            uint64_t *off = malloc((NSHARD + 1) * 8); off[0] = 0;
            for (int sh = 0; sh < NSHARD; sh++) off[sh + 1] = off[sh] + lev[pw][sh].n;
            uint64_t np = off[NSHARD];
            #pragma omp parallel
            {
                uint32_t *Bs = malloc((1u << 20) * 4);
                #pragma omp for schedule(dynamic, 16)
                for (uint64_t idx = 0; idx < np; idx++) {
                    int lo = 0, hi = NSHARD; while (hi - lo > 1) { int mid = (lo + hi) / 2; if (off[mid] <= idx) lo = mid; else hi = mid; }
                    const Key *pk = &lev[pw][lo].list[idx - off[lo]];
                    int nvp = 0; { uint32_t r0[NVMAX]; int c0[NVMAX]; unpack(pk, NVMAX, r0, c0); for (int v = 0; v < NVMAX; v++) if (r0[v]) nvp = v + 1; }
                    uint32_t up[NVMAX]; int col[NVMAX]; form_to_up(pk, nvp, up, col);
                    int order[NVMAX]; for (int l = 0; l < nvp; l++) order[l] = l;
                    for (int i = 1; i < nvp; i++) { int x = order[i], j = i; while (j > 0 && __builtin_popcount(up[order[j - 1]]) > __builtin_popcount(up[x])) { order[j] = order[j - 1]; j--; } order[j] = x; }
                    uint32_t nB = 0; upsets_rec(up, order, nvp, 0, 0, Bs, &nB);
                    int y = nvp, ny = nvp + 1, nv = nvp + 2;
                    for (uint32_t bi = 0; bi < nB; bi++) for (int q = 0; q < ncol; q++) {
                        uint32_t B = Bs[bi], rows[NVMAX], crows[NVMAX]; int c[NVMAX], ccol[NVMAX]; Key k;
                        for (int l = 0; l < nvp; l++) { rows[l] = up[l] | (1u << (l ^ 1)); c[l] = col[l]; }
                        rows[y] = 1u << ny; rows[ny] = B | (1u << y);
                        for (int l = 0; l < nvp; l++) if (B >> l & 1) rows[l ^ 1] |= 1u << y;
                        c[y] = colopt[q]; c[ny] = signed_mode ? colcode(q, s - q) : s;
                        canon(rows, c, nv, crows, ccol); pack(crows, ccol, nv, &k);
                        uint64_t h = hash_key(&k); Shard *shd = &lev[w][h & (NSHARD - 1)];
                        omp_set_lock(&shd->lock); shard_insert(shd, &k, h); omp_unset_lock(&shd->lock);
                    }
                }
                free(Bs);
            }
            free(off);
        }
        for (int sh = 0; sh < NSHARD; sh++) total += lev[w][sh].n;
        printf("U %d %llu %.1f\n", w, (unsigned long long)total, omp_get_wtime() - t0); fflush(stdout);
    }
    return 0;
}
