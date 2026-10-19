/* A123234 independent check: for every pair (lambda, mu) of cycle types of S_n except (1^n, 1^n),
 * count the n x n Latin squares L fixed by one pair (alpha, beta) of those types, i.e. with
 *     L(alpha(r), beta(c)) = L(r, c)   for all r, c,
 * directly by backtracking, without the reduction used in rc_classes.py.
 *
 * Symbol permutations commute with the row/column action and act freely, so the fixed squares fall
 * into orbits of size n! with exactly one member whose first row is 0, 1, ..., n-1.  The program
 * counts those normalized squares: N(lambda, mu) = Fix(alpha, beta) / n!.
 *
 * Output: one line per pair, "lambda;mu;N" with the parts of each cycle type space-separated.
 * burnside_check.py turns the lines into a(n).
 *
 * Build: gcc -O2 -o fixcount fixcount.c      Run: ./fixcount n [all|fpf] [skip2]
 * With "fpf", only pairs of fixed-point-free types are searched (the others are proved zero in
 * the README; running with "all" checks that too).  With "skip2", the pair (2^(n/2), 2^(n/2)) is
 * not searched (at n = 10 it has about 3.8e14 normalized fixed squares, too many to list) and is
 * reported as "skipped".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 12

static int n;
static int alpha[MAXN], beta[MAXN];
static int grid[MAXN][MAXN];
static unsigned rowm[MAXN], colm[MAXN];
static int stk[MAXN * MAXN * 2], sp;
static unsigned long long count;
static int order[MAXN]; /* rows in the order the search fills them */

/* place symbol s on the whole orbit of cell (r, c); return 0 (and undo) on a conflict */
static int place(int r, int c, int s) {
    int base = sp;
    unsigned bit = 1u << s;
    int r0 = r, c0 = c;
    do {
        if (grid[r][c] >= 0 || (rowm[r] & bit) || (colm[c] & bit)) {
            while (sp > base) {
                int cc = stk[--sp], rr = stk[--sp];
                grid[rr][cc] = -1;
                rowm[rr] &= ~bit;
                colm[cc] &= ~bit;
            }
            return 0;
        }
        grid[r][c] = s;
        rowm[r] |= bit;
        colm[c] |= bit;
        stk[sp++] = r;
        stk[sp++] = c;
        r = alpha[r];
        c = beta[c];
    } while (r != r0 || c != c0);
    return 1;
}

static void unplace(int base) {
    while (sp > base) {
        int cc = stk[--sp], rr = stk[--sp];
        unsigned bit = 1u << grid[rr][cc];
        grid[rr][cc] = -1;
        rowm[rr] &= ~bit;
        colm[cc] &= ~bit;
    }
}

static void search(int pos) {
    while (pos < n * n && grid[order[pos / n]][pos % n] >= 0) pos++;
    if (pos == n * n) {
        count++;
        return;
    }
    int r = order[pos / n], c = pos % n;
    unsigned avail = ((1u << n) - 1) & ~(rowm[r] | colm[c]);
    while (avail) {
        int s = __builtin_ctz(avail);
        avail &= avail - 1;
        int base = sp;
        if (place(r, c, s)) {
            search(pos + 1);
            unplace(base);
        }
    }
}

static void perm_of_type(const int *parts, int np, int *p) {
    int start = 0;
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < parts[i]; j++) p[start + j] = start + (j + 1) % parts[i];
        start += parts[i];
    }
}

static unsigned long long fixcount(const int *la, int na, const int *mu, int nm) {
    perm_of_type(la, na, alpha);
    perm_of_type(mu, nm, beta);
    /* Fill rows lying in short cycles of alpha first (row 0, which is normalized, stays first).
     * This only changes the search order, not the count; it makes dead pairs fail early. */
    int clen[MAXN];
    for (int r = 0; r < n; r++) {
        int k = 1;
        for (int x = alpha[r]; x != r; x = alpha[x]) k++;
        clen[r] = k;
    }
    int no = 0;
    order[no++] = 0;
    for (int k = 1; k <= n; k++)
        for (int r = 1; r < n; r++)
            if (clen[r] == k) order[no++] = r;
    memset(grid, -1, sizeof grid);
    memset(rowm, 0, sizeof rowm);
    memset(colm, 0, sizeof colm);
    sp = 0;
    count = 0;
    for (int c = 0; c < n; c++) {
        if (grid[0][c] == c) continue;
        if (grid[0][c] >= 0 || !place(0, c, c)) return 0;
    }
    search(0);
    return count;
}

/* partitions of n, parts in non-increasing order */
static int parts[64][MAXN], nparts[64], npart;
static void gen(int rem, int maxp, int *cur, int len) {
    if (rem == 0) {
        memcpy(parts[npart], cur, len * sizeof(int));
        nparts[npart++] = len;
        return;
    }
    for (int p = rem < maxp ? rem : maxp; p >= 1; p--) {
        cur[len] = p;
        gen(rem - p, p, cur, len + 1);
    }
}

static void print_type(const int *p, int k) {
    for (int i = 0; i < k; i++) printf(i ? " %d" : "%d", p[i]);
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    n = atoi(argv[1]);
    int fpf_only = argc > 2 && !strcmp(argv[2], "fpf");
    int skip2 = argc > 3 && !strcmp(argv[3], "skip2");
    int cur[MAXN];
    gen(n, n, cur, 0);
    for (int i = 0; i < npart; i++)
        for (int j = 0; j < npart; j++) {
            if (nparts[i] == n && nparts[j] == n) continue; /* identity pair: L(n)/n!, known */
            if (fpf_only && (parts[i][nparts[i] - 1] == 1 || parts[j][nparts[j] - 1] == 1)) continue;
            if (skip2 && parts[i][0] == 2 && parts[j][0] == 2 && nparts[i] == n / 2 && nparts[j] == n / 2) {
                print_type(parts[i], nparts[i]);
                printf(";");
                print_type(parts[j], nparts[j]);
                printf(";skipped\n");
                continue;
            }
            unsigned long long N = fixcount(parts[i], nparts[i], parts[j], nparts[j]);
            print_type(parts[i], nparts[i]);
            printf(";");
            print_type(parts[j], nparts[j]);
            printf(";%llu\n", N);
            fflush(stdout);
        }
    return 0;
}
