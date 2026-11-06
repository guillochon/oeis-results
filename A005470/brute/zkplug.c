/* plantri plug-in: two-sort cycle index of 3-connected planar graphs (same output as zk.py).
   Build (in the plantri source folder):  cc -o plantri_zk -O3 -DPLUGIN=zkplug.c plantri.c
   Run:  ./plantri_zk -pG n [res/mod]   (the -G switch makes plantri compute the full group)
   For every polyhedron and every automorphism g (map automorphisms, reflections included; by
   Whitney these are all graph automorphisms) it records the vertex cycle type and the edge orbits
   (length k, flipped?) of g, where flipped means g^k reverses each edge of the orbit.
   Output lines on stderr-free stdout at the end:  "<vparts>|<eorbits> <nbtot> <count>"
   meaning count automorphisms with that monomial, in graphs whose group has order nbtot, so the
   cycle-index weight of the monomial is sum count/nbtot. */

#define FILTER zk_filter
#define SUMMARY zk_summary

#define ZK_HMAX (1 << 16)
static char *zk_key[ZK_HMAX];
static long long zk_cnt[ZK_HMAX];
static int zk_grp[ZK_HMAX];
static long long zk_ident[2 * MAXE + 2];

static unsigned long zk_hash(const char *s, int g)
{
    unsigned long h = 1469598103934665603UL ^ (unsigned long)g;
    while (*s) h = (h ^ (unsigned char)*s++) * 1099511628211UL;
    return h;
}

static void zk_add(const char *key, int g)
{
    unsigned long h = zk_hash(key, g) & (ZK_HMAX - 1);
    while (zk_key[h]) {
        if (zk_grp[h] == g && strcmp(zk_key[h], key) == 0) { zk_cnt[h]++; return; }
        h = (h + 1) & (ZK_HMAX - 1);
    }
    zk_key[h] = strdup(key);
    zk_grp[h] = g;
    zk_cnt[h] = 1;
}

static int zk_cmpdesc(const void *a, const void *b) { return *(const int *)b - *(const int *)a; }
static int zk_cmpasc(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

static int zk_filter(int nbtot, int nbop, int doflip)
{
    int i, j, k, x, d;
    zk_ident[nbtot]++;
    if (nbtot == 1) return 0;
    for (j = 0; j < ne; j++) numbering[0][j]->index = j;
    int inv[MAXE], m[MAXE], seen[MAXE], vp[MAXN], vseen[MAXN], parts[MAXN], eo[MAXE];
    for (j = 0; j < ne; j++) inv[j] = numbering[0][j]->invers->index;
    for (i = 1; i < nbtot; i++) {
        for (j = 0; j < ne; j++) {
            m[j] = numbering[i][j]->index;
            vp[numbering[0][j]->start] = numbering[i][j]->start;
        }
        int np = 0;
        for (j = 0; j < nv; j++) vseen[j] = 0;
        for (j = 0; j < nv; j++) {
            if (vseen[j]) continue;
            k = 0;
            for (x = j; !vseen[x]; x = vp[x]) { vseen[x] = 1; k++; }
            parts[np++] = k;
        }
        qsort(parts, np, sizeof(int), zk_cmpdesc);
        int nedge = 0;
        for (j = 0; j < ne; j++) seen[j] = 0;
        for (d = 0; d < ne; d++) {
            if (seen[d]) continue;
            k = 0;
            x = d;
            for (;;) {
                seen[x] = 1;
                seen[inv[x]] = 1;
                x = m[x];
                k++;
                if (x == d || x == inv[d]) break;
            }
            eo[nedge++] = 2 * k + (x == inv[d]);
        }
        qsort(eo, nedge, sizeof(int), zk_cmpasc);
        char buf[4096];
        int p = 0;
        for (j = 0; j < np; j++) p += sprintf(buf + p, "%s%d", j ? "," : "", parts[j]);
        buf[p++] = '|';
        for (j = 0; j < nedge; j++) p += sprintf(buf + p, "%s%d.%d", j ? "," : "", eo[j] >> 1, eo[j] & 1);
        buf[p] = 0;
        zk_add(buf, nbtot);
    }
    return 0;
}

static void zk_summary(void)
{
    int j, g;
    for (g = 0; g <= 2 * MAXE + 1; g++)
        if (zk_ident[g]) printf("ID %d %lld\n", g, zk_ident[g]);
    for (j = 0; j < ZK_HMAX; j++)
        if (zk_key[j]) printf("%s %d %lld\n", zk_key[j], zk_grp[j], zk_cnt[j]);
    fflush(stdout);
}
