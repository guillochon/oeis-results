// For each simple matroid on k<=9 points in the rankLine catalogue, compute Aut(M) by backtracking
// and the number of Aut-orbits on compositions of s into k positive parts, for s = 1..NMAX.
// Output: totals T[s][r] = sum over simple M (k points, rank r) of orbit counts; also counts of simple M.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NMAX 10
typedef unsigned __int128 u128;
static int k, rk[512], perm[9], used[9];
static unsigned long long inv[9][10];
static long double dummy;
// group accumulation: sum over g of fixed-composition counts (exact integers)
static long long fixsum[NMAX + 1];
static long long gorder;
static void accumulate(void) {
  // cycles of perm
  int seen[9] = {0}, len[9], nc = 0;
  for (int i = 0; i < k; i++) if (!seen[i]) { int l = 0, j = i; while (!seen[j]) { seen[j] = 1; j = perm[j]; l++; } len[nc++] = l; }
  // coefficient of x^s in prod x^l/(1-x^l)
  long long poly[NMAX + 1]; memset(poly, 0, sizeof poly); poly[0] = 1;
  for (int c = 0; c < nc; c++) {
    long long np[NMAX + 1]; memset(np, 0, sizeof np);
    int l = len[c];
    for (int a = 0; a <= NMAX; a++) if (poly[a]) for (int b = a + l; b <= NMAX; b += l) np[b] += poly[a];
    memcpy(poly, np, sizeof poly);
  }
  for (int s = 0; s <= NMAX; s++) fixsum[s] += poly[s];
  gorder++;
}
static void bt(int i) {
  if (i == k) { accumulate(); return; }
  for (int t = 0; t < k; t++) {
    if (used[t]) continue;
    if (memcmp(inv[i], inv[t], sizeof inv[i])) continue;
    perm[i] = t;
    // check subsets of {0..i} containing i
    int ok = 1, base = 1 << i;
    for (int S = 0; S < base && ok; S++) {
      int X = S | base, Y = 0;
      for (int e = 0; e <= i; e++) if (X >> e & 1) Y |= 1 << perm[e];
      if (rk[X] != rk[Y]) ok = 0;
    }
    if (!ok) continue;
    used[t] = 1; bt(i + 1); used[t] = 0;
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r");
  static char line[1024];
  static long long tot[NMAX + 1][11]; static long long nsimple[10][10];
  long long lines = 0;
  while (fgets(line, sizeof line, f)) {
    lines++;
    int L = strcspn(line, "\r\n"); line[L] = 0;
    if (!strcmp(line, "NULL")) { k = 0; }
    else { k = 0; while ((1 << k) < L) k++; }
    for (int i = 0; i < (1 << k); i++) rk[i] = k ? line[i] - '0' : 0;
    int r = rk[(1 << k) - 1];
    int simple = 1;
    for (int e = 0; e < k && simple; e++) if (rk[1 << e] != 1) simple = 0;
    for (int a = 0; a < k && simple; a++) for (int b = a + 1; b < k; b++) if (rk[(1 << a) | (1 << b)] != 2) { simple = 0; break; }
    if (!simple) continue;
    nsimple[k][r]++;
    // invariants: for each element, number of subsets containing it of each rank, weighted by size
    memset(inv, 0, sizeof inv);
    for (int X = 0; X < (1 << k); X++) { int sz = __builtin_popcount(X); for (int e = 0; e < k; e++) if (X >> e & 1) inv[e][rk[X]] += 1ULL << (5 * sz % 60); }
    memset(fixsum, 0, sizeof fixsum); gorder = 0; memset(used, 0, sizeof used);
    bt(0);
    for (int s = 0; s <= NMAX; s++) {
      if (fixsum[s] % gorder) { fprintf(stderr, "nonint line %lld\n", lines); return 1; }
      tot[s][r] += fixsum[s] / gorder;
    }
  }
  fprintf(stderr, "lines %lld\n", lines);
  for (int kk = 0; kk <= 9; kk++) { printf("simple k=%d:", kk); for (int r = 0; r <= 9; r++) printf(" %lld", nsimple[kk][r]); printf("\n"); }
  for (int s = 0; s <= NMAX; s++) { printf("loopless-from-k<=9 s=%d:", s); for (int r = 0; r <= 10; r++) printf(" %lld", tot[s][r]); printf("\n"); }
  return 0;
}
