// Compact version of minregion_tm.cpp (same method; see the comment there).  States are packed into 192 bits:
// 4-bit component labels for up to 20 frontier cells, then per component 9 bits of capped orbit-size sum S and
// 2 wall-flag bits, up to 10 components.  Counts are kept both mod 2^128 and mod p = 2^61 - 1, so T(n,k) is
// recovered exactly (CRT) up to 2^189, enough for octants of up to 189 cells (n <= 37).
// Usage: minregion_tm2 n k   prints "n k T(n,k)".
#include <bits/stdc++.h>
using namespace std;
typedef unsigned __int128 u128;
typedef unsigned long long u64;
static const u64 P61 = (1ULL << 61) - 1;
static inline u64 addp(u64 a, u64 b) { u64 s = a + b; if (s >= P61) s -= P61; return s; }

struct Key { u64 w[3]; bool operator==(const Key &o) const { return w[0] == o.w[0] && w[1] == o.w[1] && w[2] == o.w[2]; } };
struct Ent { Key k; u128 c; u64 m; };
static inline u64 hk(const Key &k) { u64 x = k.w[0] * 0x9E3779B97F4A7C15ULL ^ (k.w[1] + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL ^ k.w[2] * 0x165667B19E3779F9ULL; x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ULL; return x ^ (x >> 29); }

// Open-addressing table in one malloc block; a stored key has bit 63 of w[2] set (bits 190-191 are unused), so
// an all-zero slot is empty.  Grows by 1.5x at load 0.85; after a step it is compacted in place and shrunk with
// realloc, becoming the next step's state list (no second copy).
static const u64 USED = 1ULL << 63;
struct Table {
    Ent *e = nullptr; u64 cap = 0, cnt = 0;
    static Ent *alloc(u64 c) { Ent *p = (Ent *)calloc(c, sizeof(Ent)); if (!p) { fprintf(stderr, "out of memory\n"); exit(2); } return p; }
    void init(u64 want) { cap = max<u64>(1024, (u64)(want / 0.8)); e = alloc(cap); cnt = 0; }
    inline u64 slot(const Key &k) const { return (u64)(((u128)hk(k) * cap) >> 64); }
    void insert_raw(const Ent &x) { u64 i = slot(x.k); while (e[i].k.w[2]) { if (++i == cap) i = 0; } e[i] = x; }
    void grow() { Ent *old = e; u64 oc = cap; cap = cap + cap / 2; e = alloc(cap); for (u64 i = 0; i < oc; i++) if (old[i].k.w[2]) insert_raw(old[i]); free(old); }
    void put(Key k, u128 c, u64 m) {
        k.w[2] |= USED;
        if (cnt + 1 > cap * 0.85) grow();
        u64 i = slot(k);
        while (e[i].k.w[2]) { if (e[i].k == k) { e[i].c += c; e[i].m = addp(e[i].m, m); return; } if (++i == cap) i = 0; }
        e[i].k = k; e[i].c = c; e[i].m = m; cnt++;
    }
    Ent *compact(u64 &n) { u64 j = 0; for (u64 i = 0; i < cap; i++) if (e[i].k.w[2]) e[j++] = e[i]; Ent *p = (Ent *)realloc(e, max<u64>(j, 1) * sizeof(Ent)); n = j; e = nullptr; cap = cnt = 0; return p ? p : e; }
};

static string str128(u128 x) { if (!x) return "0"; string s; while (x) { s += char('0' + int(x % 10)); x /= 10; } reverse(s.begin(), s.end()); return s; }

int n, K, h, NV;
static inline int Hsize(int f) { return f == 0 ? 1 : (f == 3 ? 8 : 2); }
static inline int capFor(int f) { int hs = Hsize(f); return (8 * K + hs - 1) / hs; }
static inline int getb(const Key &k, int pos, int len) { int wi = pos >> 6, of = pos & 63; u64 v = k.w[wi] >> of; if (of + len > 64) v |= k.w[wi + 1] << (64 - of); return int(v & ((1ULL << len) - 1)); }
static inline void setb(Key &k, int pos, int len, u64 v) { int wi = pos >> 6, of = pos & 63; k.w[wi] |= v << of; if (of + len > 64) k.w[wi + 1] |= v >> (64 - of); }

int main(int argc, char **argv) {
    n = atoi(argv[1]); K = atoi(argv[2]); h = (n + 1) / 2;
    if (8 * K > 511) { fprintf(stderr, "k too large for 9-bit sums\n"); return 1; }
    map<pair<int,int>, int> id; vector<pair<int,int>> cell;
    for (int j = 0; j < h; j++) for (int i = 0; i <= j; i++) { id[{i, j}] = cell.size(); cell.push_back({i, j}); }
    NV = cell.size();
    if (NV > 189) { fprintf(stderr, "octant too large for exact CRT\n"); return 1; }
    auto rep = [&](int r, int c) { int a = min(r, n - 1 - r), b = min(c, n - 1 - c); return id[{min(a, b), max(a, b)}]; };
    vector<int> orb_s(NV, 0), fD(NV), fV(NV), lastNb(NV);
    for (int r = 0; r < n; r++) for (int c = 0; c < n; c++) orb_s[rep(r, c)]++;
    vector<set<int>> adj(NV);
    auto near = [&](int r1, int c1, int r2, int c2) { return abs(r1 - r2) <= 1 && abs(c1 - c2) <= 1; };
    for (int v = 0; v < NV; v++) {
        int r = cell[v].first, c = cell[v].second;
        fD[v] = near(r, c, c, r); fV[v] = near(r, c, r, n - 1 - c);
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
            if (!a && !b) continue; int x = r + a, y = c + b;
            if (x < 0 || y < 0 || x >= n || y >= n) continue;
            int w = rep(x, y); if (w != v) { adj[v].insert(w); adj[w].insert(v); }
        }
    }
    vector<vector<int>> nbBefore(NV);
    for (int v = 0; v < NV; v++) { lastNb[v] = v; for (int w : adj[v]) { lastNb[v] = max(lastNb[v], w); if (w < v) nbBefore[v].push_back(w); } }

    vector<int> front;
    u64 ncur = 1; Ent *cur = (Ent *)calloc(1, sizeof(Ent)); cur[0].c = 1; cur[0].m = 1;
    Table nxt; size_t peak = 0, prevsz = 1; double growth = 2.0;
    for (int t = 0; t < NV; t++) {
        vector<int> pos; for (int w : nbBefore[t]) { int p = find(front.begin(), front.end(), w) - front.begin(); if (p == (int)front.size()) { fprintf(stderr, "frontier bug\n"); return 1; } pos.push_back(p); }
        vector<int> nf = front; nf.push_back(t);
        vector<int> keep; for (int p = 0; p < (int)nf.size(); p++) if (lastNb[nf[p]] > t) keep.push_back(p);
        int F = front.size();
        if ((int)keep.size() > 20) { fprintf(stderr, "frontier too wide\n"); return 1; }
        nxt.init((u64)(ncur * growth));
        int sbase = 4 * F;   // component data starts after the labels of the current frontier
        // consume the state list from the end, giving memory back as it shrinks
        for (u64 ii = ncur; ii-- > 0;) {
            if (ii > 0 && ii % (1u << 22) == 0) { Ent *q = (Ent *)realloc(cur, (ii + 1) * sizeof(Ent)); if (q) cur = q; }
            const Ent E = cur[ii];
            int lab[24], nc = 0, S[24], fl[24];
            for (int p = 0; p < F; p++) { lab[p] = getb(E.k, 4 * p, 4); nc = max(nc, lab[p]); }
            for (int c = 1; c <= nc; c++) { int v = getb(E.k, sbase + 11 * (c - 1), 11); S[c] = v & 511; fl[c] = v >> 9; }
            for (int x = 0; x <= 1; x++) {
                int L[24], SS[24], FF[24], ncc = nc;
                for (int p = 0; p < F; p++) L[p] = lab[p];
                for (int c = 1; c <= nc; c++) { SS[c] = S[c]; FF[c] = fl[c]; }
                if (!x) L[F] = 0;
                else {
                    int nw = ++ncc; SS[nw] = orb_s[t]; FF[nw] = fD[t] | (fV[t] << 1);
                    bool merged[24] = {false};
                    for (int p : pos) { int c = L[p]; if (c && c != nw && !merged[c]) { merged[c] = true; SS[nw] += SS[c]; FF[nw] |= FF[c]; } }
                    for (int p = 0; p < F; p++) if (L[p] && merged[L[p]]) L[p] = nw;
                    L[F] = nw; SS[nw] = min(SS[nw], capFor(FF[nw]));
                }
                char alive[24] = {0};
                for (int p : keep) if (L[p]) alive[L[p]] = 1;
                bool ok = true;
                for (int p = 0; p <= F && ok; p++) { int c = L[p]; if (!c || alive[c]) continue; if (SS[c] * Hsize(FF[c]) < 8 * K) ok = false; alive[c] = 2; }
                if (!ok) continue;
                Key nk{{0, 0, 0}}; int mp[24] = {0}, m = 0, order[24];
                for (size_t q = 0; q < keep.size(); q++) { int c = L[keep[q]]; if (c && !mp[c]) { mp[c] = ++m; order[m - 1] = c; } setb(nk, 4 * q, 4, c ? mp[c] : 0); }
                if (m > 10 || m > 15) { fprintf(stderr, "too many components\n"); exit(1); }
                int nb = 4 * keep.size();
                for (int q = 0; q < m; q++) { int c = order[q]; int v = min(SS[c], capFor(FF[c])); setb(nk, nb + 11 * q, 11, (u64)(v | (FF[c] << 9))); }
                nxt.put(nk, E.c, E.m);
            }
        }
        free(cur); u64 prev = ncur; cur = nxt.compact(ncur); (void)prev;
        vector<int> f2; for (int p : keep) f2.push_back(nf[p]); front = f2;
        growth = min(2.0, max(1.0, (double)ncur / max<size_t>(1, prevsz)) * 1.05); prevsz = ncur;
        peak = max<size_t>(peak, ncur);
    }
    u128 c = 0; u64 m = 0; for (u64 i = 0; i < ncur; i++) { c += cur[i].c; m = addp(m, cur[i].m); }
    // CRT: T = c + 2^128 * t with T = m (mod P61); 2^128 mod P61 = 2^(128 mod 61) = 2^6
    u64 cm = (u64)(c % P61); u64 inv64 = 0; { u64 a = 64, e = P61 - 2, r = 1; while (e) { if (e & 1) r = (u128)r * a % P61; a = (u128)a * a % P61; e >>= 1; } inv64 = r; }
    u64 tt = (u128)((m + P61 - cm) % P61) * inv64 % P61;
    string out;
    if (tt == 0) out = str128(c);
    else { // print c + tt * 2^128 in decimal via base-1e18 limbs
        vector<u64> d; // little-endian base 1e18
        auto mul_add = [&](vector<u64> &a, u64 mul, u64 add) { u128 carry = add; for (auto &x : a) { u128 v = (u128)x * mul + carry; x = (u64)(v % 1000000000000000000ULL); carry = v / 1000000000000000000ULL; } while (carry) { a.push_back((u64)(carry % 1000000000000000000ULL)); carry /= 1000000000000000000ULL; } };
        d.push_back(0); mul_add(d, 1, tt);
        for (int i = 0; i < 4; i++) mul_add(d, 1ULL << 32, 0);
        u128 cc = c; vector<u64> cl; while (cc) { cl.push_back((u64)(cc % 1000000000000000000ULL)); cc /= 1000000000000000000ULL; }
        u128 carry = 0; for (size_t i = 0; i < max(d.size(), cl.size()) || carry; i++) { if (i >= d.size()) d.push_back(0); u128 v = (u128)d[i] + (i < cl.size() ? cl[i] : 0) + carry; d[i] = (u64)(v % 1000000000000000000ULL); carry = v / 1000000000000000000ULL; }
        char buf[32]; snprintf(buf, sizeof buf, "%llu", (unsigned long long)d.back()); out = buf;
        for (int i = (int)d.size() - 2; i >= 0; i--) { snprintf(buf, sizeof buf, "%018llu", (unsigned long long)d[i]); out += buf; }
    }
    printf("%d %d %s\n", n, K, out.c_str());
    fprintf(stderr, "n=%d k=%d octant=%d peak_states=%zu\n", n, K, NV, peak);
    return 0;
}
