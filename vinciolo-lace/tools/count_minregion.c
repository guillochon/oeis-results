/* Count n x n 0/1 arrays (mode 0: D4-symmetric, mode 1: all) in which every 8-connected region of 1s has
   at least k cells, for k = 1..8 at once.  (k = 2 is "no pinhole" in the lace reading: 1 = open mesh.) */
#include <stdio.h>
#include <stdlib.h>
int n, A[16][16], seen[16][16], qi[256], qj[256];
int main(int argc, char **argv) {
  n = atoi(argv[1]); int mode = atoi(argv[2]); int id[16][16], K = 0, oid[8][8];
  for (int a = 0; a < 8; a++) for (int b = 0; b < 8; b++) oid[a][b] = -1;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
    if (mode) { id[i][j] = K++; continue; }
    int x = i < n-1-i ? i : n-1-i, y = j < n-1-j ? j : n-1-j, a = x < y ? x : y, b = x < y ? y : x;
    if (oid[a][b] < 0) oid[a][b] = K++; id[i][j] = oid[a][b]; }
  unsigned long long cnt[10] = {0};
  for (unsigned long long m = 0; m < (1ULL << K); m++) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { A[i][j] = (m >> id[i][j]) & 1; seen[i][j] = 0; }
    int mn = 1000;
    for (int i = 0; i < n && mn > 1; i++) for (int j = 0; j < n; j++) if (A[i][j] && !seen[i][j]) {
      int h = 0, t = 0; qi[t] = i; qj[t++] = j; seen[i][j] = 1;
      while (h < t) { int u = qi[h], v = qj[h++];
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) { int x = u + a, y = v + b;
          if (x >= 0 && y >= 0 && x < n && y < n && A[x][y] && !seen[x][y]) { seen[x][y] = 1; qi[t] = x; qj[t++] = y; } } }
      if (t < mn) mn = t; if (mn == 1) break; }
    for (int k = 1; k <= 8; k++) if (mn >= k) cnt[k]++; }
  printf("n=%d mode=%d", n, mode); for (int k = 1; k <= 8; k++) printf(" %llu", cnt[k]); printf("\n");
}
