// Count n x n 0/1 designs symmetric under the dihedral group of the square (D4) with no pinhole:
// no 0 (open) cell whose 8 neighbors are all 1 (filled); cells outside the grid count as filled (the frame).
//
// Method: a D4-symmetric design is fixed by its values on the octant cells (a, b), 0 <= a <= b < h,
// h = ceil(n/2). Two octant cells are neighbors when some full-grid cells in their orbits are king-adjacent;
// the pinhole condition for every full cell in an orbit is the same as for one representative, and an open
// representative is satisfied iff some neighbor of it (in the full grid) is open. We place octant cells in
// order (b ascending, then a) and keep a frontier: each placed cell still needed later stores
// (value, satisfied). A cell leaves the frontier once all its neighbors are placed; it must then be filled
// or satisfied. Counts are exact multi-limb integers (additions only).
//
//   pinhole_tm n              prints n, a(n), max frontier states, frontier slots
//   pinhole_tm_mod n P        same, a(n) mod P (P = 0: mod 2^64); compiled with -DMODP
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#ifdef MODP
// residue mode: counts mod P (P = 0 means mod 2^64, plain wrapping); recombine runs with CRT.
static uint64_t P = 0;
struct Big { uint64_t w[1]; };
static inline void add(Big &a, const Big &b) {
    if (!P) { a.w[0] += b.w[0]; return; }
    uint64_t s = a.w[0] + b.w[0]; if (s >= P) s -= P; a.w[0] = s;     // P < 2^63, no overflow
}
static const int L = 1;
#else
static const int L = 6;                       // 384-bit counts: enough for n <= 47 (T(24) = 300 bits)
struct Big { uint64_t w[L]; };
static inline void add(Big &a, const Big &b) {
    unsigned __int128 c = 0;
    for (int i = 0; i < L; i++) { c += (unsigned __int128)a.w[i] + b.w[i]; a.w[i] = (uint64_t)c; c >>= 64; }
}
#endif
static string dec(Big a) {
    string s;
    while (true) {
        bool zero = true; for (int i = 0; i < L; i++) if (a.w[i]) zero = false;
        if (zero) break;
        unsigned __int128 r = 0;
        for (int i = L - 1; i >= 0; i--) { unsigned __int128 cur = (r << 64) | a.w[i]; a.w[i] = (uint64_t)(cur / 10); r = cur % 10; }
        s += char('0' + (int)r);
    }
    if (s.empty()) s = "0";
    reverse(s.begin(), s.end()); return s;
}

// open-addressing hash map: key (uint64) -> Big
struct Map {
    vector<uint64_t> keys; vector<Big> vals; vector<uint8_t> used; size_t cnt = 0, mask = 0;
    void init(size_t cap) { size_t m = 16; while (m < 2 * cap) m <<= 1; keys.assign(m, 0); vals.assign(m, Big{}); used.assign(m, 0); mask = m - 1; cnt = 0; }
    void grow() {
        Map nm; nm.init(cnt * 2 + 16);
        for (size_t i = 0; i <= mask; i++) if (used[i]) nm.put(keys[i], vals[i]);
        *this = std::move(nm);
    }
    void put(uint64_t k, const Big &v) {
        if (2 * (cnt + 1) > mask + 1) grow();
        size_t i = (k * 0x9E3779B97F4A7C15ULL) >> 20 & mask;
        while (used[i] && keys[i] != k) i = (i + 1) & mask;
        if (!used[i]) { used[i] = 1; keys[i] = k; vals[i] = v; cnt++; } else add(vals[i], v);
    }
};

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
#ifdef MODP
    if (argc > 2) P = strtoull(argv[2], 0, 10);
#endif
    int h = (n + 1) / 2;
    // octant cells in placement order
    vector<pair<int,int>> cells; vector<vector<int>> id(h, vector<int>(h, -1));
    for (int b = 0; b < h; b++) for (int a = 0; a <= b; a++) { id[a][b] = cells.size(); cells.push_back({a, b}); }
    int K = cells.size();
    auto fold = [&](int x) { return min(x, n - 1 - x); };          // full index -> quarter index
    // neighbors of each octant cell, via its representative full cell (a, b)
    vector<vector<int>> nb(K);
    for (int c = 0; c < K; c++) {
        auto [a, b] = cells[c];
        for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++) {
            if (!di && !dj) continue;
            int i = a + di, j = b + dj;
            if (i < 0 || j < 0 || i >= n || j >= n) continue;          // frame: filled, never helps
            int x = fold(i), y = fold(j); if (x > y) swap(x, y);
            int d = id[x][y];
            if (find(nb[c].begin(), nb[c].end(), d) == nb[c].end()) nb[c].push_back(d);
        }
    }
    // last time each cell is needed
    vector<int> last(K);
    for (int c = 0; c < K; c++) { last[c] = c; for (int d : nb[c]) last[c] = max(last[c], d); }
    // frontier slot assignment: slots are reused; state packs 2 bits per slot (bit0 = open, bit1 = satisfied)
    vector<int> slot(K, -1); vector<int> freeSlots; int nslots = 0;
    vector<vector<int>> retire(K);                                   // cells retiring right after placing c
    for (int c = 0; c < K; c++) retire[last[c]].push_back(c);
    vector<int> slotOf(K);
    {
        vector<int> fr; int used = 0;
        for (int c = 0; c < K; c++) {
            int s;
            if (!freeSlots.empty()) { s = freeSlots.back(); freeSlots.pop_back(); } else s = used++;
            slotOf[c] = s;
            for (int d : retire[c]) freeSlots.push_back(slotOf[d]);
        }
        nslots = used;
    }
    if (2 * nslots > 64) { fprintf(stderr, "frontier too wide: %d slots\n", nslots); return 1; }

    Map cur; cur.init(1024); Big one{}; one.w[0] = 1; cur.put(0, one);
    size_t maxStates = 1;
    for (int c = 0; c < K; c++) {
        Map nxt; nxt.init(cur.cnt + 16);
        int sc = slotOf[c];
        bool selfNb = find(nb[c].begin(), nb[c].end(), c) != nb[c].end();
        vector<int> earlier; for (int d : nb[c]) if (d < c) earlier.push_back(d);
        for (size_t i = 0; i <= cur.mask; i++) {
            if (!cur.used[i]) continue;
            uint64_t st = cur.keys[i];
            for (int v = 0; v <= 1; v++) {                           // v = 1: open
                uint64_t s = st & ~(3ULL << (2 * sc));
                if (v) {
                    bool sat = selfNb;
                    for (int d : earlier) {
                        int sd = slotOf[d];
                        if ((s >> (2 * sd)) & 1) { sat = true; s |= 2ULL << (2 * sd); }   // both see each other
                    }
                    s |= (1ULL | (sat ? 2ULL : 0)) << (2 * sc);
                }
                bool ok = true;
                for (int d : retire[c]) {
                    int sd = slotOf[d]; uint64_t f = (s >> (2 * sd)) & 3;
                    if (f == 1) { ok = false; break; }                // open, never satisfied: pinhole
                    s &= ~(3ULL << (2 * sd));
                }
                if (ok) nxt.put(s, cur.vals[i]);
            }
        }
        cur = std::move(nxt);
        maxStates = max(maxStates, cur.cnt);
    }
    Big tot{};
    for (size_t i = 0; i <= cur.mask; i++) if (cur.used[i]) add(tot, cur.vals[i]);
    printf("%d %s %zu %d\n", n, dec(tot).c_str(), maxStates, nslots);
    return 0;
}
