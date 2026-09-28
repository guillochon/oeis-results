// a(8) of OEIS A027624 (independent sets of Q_8) on the GPU.
//
// Union grouping (see ../rust/src/ygroup.rs, which this mirrors term by term):
//
//     a(8) = sum over Y subset of V(Q_6) of c(Y) * f_6(~Y)^2,
//     c(Y) = 3^(isolated vertices of Q_6[Y]) * 2^(other components of Q_6[Y]).
//
// Q_6 = Q_4 x C_4, four 16-bit fibres Y = Y0 | Y1<<16 | Y2<<32 | Y3<<48. Y_lo = Y0 | Y1<<16 runs
// over the 1,228,158 orbit representatives of Aut(Q_5) (file from `hypercube-indep reps6`),
// Y3 and Y2 over all 2^16 values. One block = one (rep, Y3) slice = 2^16 terms:
//
//     f_6(W) = sum_{t2 in I_4, t2 inside W2} G(t2),
//     G(t2)  = sum_{t0 in I_4, t0 inside W0} f_4(W1 \ (t0|t2)) * f_4(W3 \ (t0|t2)).
//
// Every I_4 set is l | h<<8 with l, h in I_3 (35 sets), so G is stored as a 35 x 35 matrix. Each warp
// takes one high byte H of W2 at a time (handed out dynamically, most work first):
// g[l] = sum over h inside H of G[l][h] (the h loop is the same for every lane), then an 8-bit
// subset-sum transform in registers (warp shuffles) gives f_6 for the 256 low bytes.
// Terms are accumulated exactly in 192 bits (every term is at most a(8) < 2^131).
// Tuned with Nsight Compute (see ../README.md): < 33 KB shared memory so 3 blocks fit per SM.
//
// Build (WSL):  /usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -o a8 a8.cu
// Usage:        ./a8 REPS.bin --from I --to J [--stride S] [--batch B] [--out FILE] [--sym]
//               ./a8 REPS.bin --check REFERENCE.txt [--sym]   (compare with `hypercube-indep cpu6 [--sym]`)
// --sym uses the lo <-> hi swap symmetry (about 45% fewer terms); the per-rep sums are then the
// weighted ones, and must be compared with `cpu6 --sym`.
// Output lines: rep_index rep_mask orbit_size unweighted_sum
// Rerunning with the same --out resumes: indices already in the file are skipped.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <chrono>
#include <thread>
#include <unistd.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(1); } } while (0)

static const int NI4 = 743;          // independent sets of Q_4
static const int NI3 = 35;           // independent sets of Q_3 (possible bytes of an I_4 set)
#ifndef THREADS
#define THREADS 256                  // 8 warps; shared memory is kept under 33 KB so 3 blocks fit per SM
#endif
#define NWARP (THREADS / 32)


__constant__ u64 c_p3lo[65], c_p3hi[65]; // 3^i as 128-bit
__constant__ u16 c_n4[16];               // Q_4 neighbours of each fibre position
__constant__ u64 c_hmask[256];           // bit h set iff I_3[h] is inside the byte H
__constant__ u8 c_i3[NI3];               // the independent sets of Q_3 (as bytes)
__constant__ u8 c_i3idx[256];            // index of a byte in I_3, or 255
__constant__ u8 c_border[256];           // bytes sorted by popcount (so skipped terms come in whole warps)
__constant__ u8 c_aorder[256];           // bytes by decreasing popcount (most work first, for balance)

// ------------------------------------------------------------------------------------------------
// 192-bit helpers
// ------------------------------------------------------------------------------------------------
struct U192 { u64 w[3]; };

__host__ __device__ inline void add192(U192 &a, u64 p0, u64 p1, u64 p2) {
    u64 s0 = a.w[0] + p0; u64 c0 = s0 < p0;
    u64 s1 = a.w[1] + p1; u64 c1 = s1 < p1;
    u64 s1b = s1 + c0;   c1 += s1b < s1;
    a.w[0] = s0; a.w[1] = s1b; a.w[2] += p2 + c1;
}

// ------------------------------------------------------------------------------------------------
// Device code
// ------------------------------------------------------------------------------------------------
__device__ inline u64 nbrs64(u64 x) {
    u64 r = ((x & 0x5555555555555555ULL) << 1) | ((x >> 1) & 0x5555555555555555ULL);
    r |= ((x & 0x3333333333333333ULL) << 2) | ((x >> 2) & 0x3333333333333333ULL);
    r |= ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    r |= ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x >> 8) & 0x00FF00FF00FF00FFULL);
    r |= (x << 16) | (x >> 48) | (x >> 16) | (x << 48);
    return r;
}

__device__ inline u32 nbrs16(u32 y) {
    u32 r = ((y & 0x5555) << 1) | ((y >> 1) & 0x5555);
    r |= ((y & 0x3333) << 2) | ((y >> 2) & 0x3333);
    r |= ((y & 0x0F0F) << 4) | ((y >> 4) & 0x0F0F);
    r |= ((y & 0x00FF) << 8) | ((y >> 8) & 0x00FF);
    return r & 0xFFFF;
}

extern "C" __global__ void __launch_bounds__(THREADS)
slice_kernel(const u32 *__restrict__ reps, int rep0,
             const u16 *__restrict__ fm,        // f_4 table, 2^16
             const u16 *__restrict__ im,        // I_4 (743)
             const u8 *__restrict__ hid,        // index in I_3 of the high byte of each I_4 set
             const u8 *__restrict__ lid,        // index in I_3 of the low byte of each I_4 set
             const u32 *__restrict__ lstoff,    // CSR offsets over X (2^16 + 1)
             const u16 *__restrict__ lst,       // I_4 sets inside X
             u64 *__restrict__ out,             // 3 limbs per block
             int sym)                           // use the lo <-> hi swap symmetry (see below)
{
    // GdT[h][l] = G(I_3[l] | I_3[h] << 8) (zero when that is not an I_4 set), stored transposed so
    // that lanes (indexed by l) read consecutive words.
    __shared__ u64 GdT[NI3 * NI3];
    __shared__ u64 gz[NWARP][256];     // per warp: f_6 for the 256 low bytes (first 35 words reused
                                       // for the filtered sums before the transform)
    __shared__ u32 nb[2][256];         // extended neighbours of fibre-2 positions, by byte
    __shared__ u32 TL[256], TH[256];   // bit q: fixed component q does not touch the low / high byte
    __shared__ u64 p3lo[65], p3hi[65]; // lane-indexed tables in shared memory (constant memory
    __shared__ u8 border[256], i3idx[256], aorder[256]; // serialises divergent reads)
    __shared__ u8 i3s[NI3];
    __shared__ u32 A[16], touch2[32];
    __shared__ u32 s_nt, s_y13, s_i1, s_i3, s_base_iso, s_base_c2, s_next;
    __shared__ u64 red[3][NWARP];

    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const u32 ylo = reps[rep0 + blockIdx.y];
    const u32 y0 = ylo & 0xFFFF, y1 = ylo >> 16, y3 = blockIdx.x;
    const u32 w0 = ~y0 & 0xFFFF, w1 = ~y1 & 0xFFFF, w3 = ~y3 & 0xFFFF;

    // Swap symmetry: the summand is invariant under exchanging Y_lo and Y_hi = (Y3, Y2), and
    // popcount is invariant under Aut(Q_5). So terms with |Y_hi| < |Y_lo| are skipped, equal ones
    // count once and larger ones twice (the weighted total is unchanged).
    const u32 pl = __popc(ylo), p3 = __popc(y3);
    if (sym && p3 + 16 < pl) {
        if (tid == 0) { u64 *o = out + 3 * ((u64)blockIdx.y * gridDim.x + blockIdx.x); o[0] = o[1] = o[2] = 0; }
        return;
    }

    for (int e = tid; e < NI3 * NI3; e += THREADS) GdT[e] = 0;
    for (int e = tid; e < 256; e += THREADS) { border[e] = c_border[e]; i3idx[e] = c_i3idx[e]; aorder[e] = c_aorder[e]; }
    if (tid < 65) { p3lo[tid] = c_p3lo[tid]; p3hi[tid] = c_p3hi[tid]; }
    if (tid < NI3) i3s[tid] = c_i3[tid];
    if (tid == 0) s_next = 0;
    __syncthreads();

    // G(t2) over I_4
    const u32 a0 = lstoff[w0];
#ifdef SKIP_G
    const u32 a1 = a0;   // timing experiment: no G work
#else
    const u32 a1 = lstoff[w0 + 1];
#endif
    for (int idx = tid; idx < NI4; idx += THREADS) {
        u32 t2 = im[idx];
        u64 g = 0;
        for (u32 q = a0; q < a1; q++) {
            u32 u = lst[q] | t2;
            g += (u64)fm[w1 & ~u] * (u64)fm[w3 & ~u];
        }
        GdT[hid[idx] * NI3 + lid[idx]] = g;
    }

    // components of the fixed part (fibres 0, 1, 3): one thread
    if (tid == 0) {
        u64 fixed = (u64)y0 | ((u64)y1 << 16) | ((u64)y3 << 48);
        u64 iso_f = fixed & ~nbrs64(fixed);
        s_y13 = (y1 | y3) & 0xFFFF;
        s_i1 = (u32)(iso_f >> 16) & 0xFFFF;
        s_i3 = (u32)(iso_f >> 48) & 0xFFFF;
        s_base_iso = __popcll(iso_f & 0xFFFF);
        for (int p = 0; p < 16; p++) A[p] = c_n4[p];
        u32 nt = 0, bc2 = 0;
        u64 rem = fixed & ~iso_f;
        while (rem) {
            u64 comp = rem & (~rem + 1);
            while (true) {
                u64 nxt = (comp | nbrs64(comp)) & fixed;
                if (nxt == comp) break;
                comp = nxt;
            }
            rem &= ~comp;
            u32 tau = (u32)((comp >> 16) | (comp >> 48)) & 0xFFFF;
            if (tau == 0) { bc2++; }
            else {
                touch2[nt++] = tau;
                for (int p = 0; p < 16; p++) if (tau >> p & 1) A[p] |= tau;
            }
        }
        s_nt = nt; s_base_c2 = bc2;
    }
    __syncthreads();
    for (int e = tid; e < 512; e += THREADS) {
        int ch = e >> 8, b = e & 255;
        u32 v = 0;
        for (int j = 0; j < 8; j++) if (b >> j & 1) v |= A[8 * ch + j];
        nb[ch][b] = v;
        u32 t = 0;
        for (u32 q = 0; q < s_nt; q++) t |= (u32)((((touch2[q] >> (8 * ch)) & 255) & b) == 0) << q;
        (ch ? TH : TL)[b] = t;
    }
    __syncthreads();

    // the block-constant factor 3^base_iso * 2^base_c2 is applied once at the end
    const u32 y13 = s_y13, i1 = s_i1, i3 = s_i3;
    U192 acc = {{0, 0, 0}};
    u64 *g = gz[warp];

    while (true) {
        // high bytes of Y2 are handed out dynamically, most work first
        u32 k = 0;
        if (lane == 0) k = atomicAdd(&s_next, 1);
        k = __shfl_sync(0xFFFFFFFF, k, 0);
        if (k >= 256) break;
        const u32 a = aorder[k];
        const u32 p3a = p3 + __popc(a);
        if (sym && p3a + 8 < pl) break;                   // this and every later a: all terms skipped
        const u64 hm = c_hmask[~a & 255];                 // I_3 sets inside the high byte of W2 (warp-uniform)

        // the 35 filtered sums g[l] = sum over h inside H of GdT[h][l]: lane l (and lane l-32 for
        // l = 32..34, in the same loop)
        u64 s = 0, s2 = 0;
        for (u32 m = (u32)hm; m; m &= m - 1) {
            const u64 *row = GdT + (__ffs(m) - 1) * NI3;
            s += row[lane];
            if (lane < NI3 - 32) s2 += row[32 + lane];
        }
        for (u32 m = (u32)(hm >> 32); m; m &= m - 1) {
            const u64 *row = GdT + (31 + __ffs(m)) * NI3;
            s += row[lane];
            if (lane < NI3 - 32) s2 += row[32 + lane];
        }
        g[lane] = s;
        if (lane < NI3 - 32) g[32 + lane] = s2;
        __syncwarp();
        // 8-bit subset-sum transform in registers: entry l = lane + 32 j is v[j] of this lane;
        // bits 0-4 of l are lane bits (shuffles), bits 5-7 are j bits (adds)
        u64 v[8];
#pragma unroll
        for (int j = 0; j < 8; j++) { u32 li = i3idx[lane + 32 * j]; v[j] = li != 255 ? g[li] : 0; }
        __syncwarp();
#pragma unroll
        for (int bit = 1; bit < 32; bit <<= 1) {
#pragma unroll
            for (int j = 0; j < 8; j++) {
                u64 t = __shfl_xor_sync(0xFFFFFFFF, v[j], bit);
                if (lane & bit) v[j] += t;
            }
        }
#pragma unroll
        for (int jb = 1; jb < 8; jb <<= 1) {
#pragma unroll
            for (int j = 0; j < 8; j++) if (j & jb) v[j] += v[j ^ jb];
        }
#pragma unroll
        for (int j = 0; j < 8; j++) g[lane + 32 * j] = v[j];
        __syncwarp();

        for (int j = 0; j < 8; j++) {
#ifdef NATURAL_B
            const u32 b = lane + 32 * j;     // conflict-free g / TL reads, but skips are per lane
#else
            const u32 b = border[lane + 32 * j];
#endif
            u32 extra = 0;
            if (sym) {
                const u32 ph = p3a + __popc(b);
                if (ph < pl) continue;
                extra = ph > pl;
            }
            const u32 y2 = (a << 8) | b;
            const u64 f = g[~b & 255];
            // c(Y): isolated vertices and components of size >= 2
            u32 isoq = y2 & ~nbrs16(y2) & ~y13;
            u32 iso = __popc(i1 & ~y2) + __popc(i3 & ~y2) + __popc(isoq);   // beyond base_iso
            u32 n = 0;                                                      // beyond base_c2
#ifdef SKIP_C
            u32 rem = 0;   // timing experiment: no component count
#else
            u32 rem = y2 & ~isoq;
#endif
            while (rem) {
                u32 comp = rem & (~rem + 1);
                while (true) {
                    u32 nxt = (comp | nb[0][comp & 255] | nb[1][comp >> 8]) & y2;
                    if (nxt == comp) break;
                    comp = nxt;
                }
                rem &= ~comp;
                n++;
            }
            n += __popc(TL[b] & TH[a]) + extra;   // untouched fixed components; weight 2
            // term = 3^iso * (f^2 << n), exact
            u64 f2lo = f * f, f2hi = __umul64hi(f, f);
            u64 x0 = f2lo << n;
            u64 x1 = n ? ((f2hi << n) | (f2lo >> (64 - n))) : f2hi;
            u64 v0 = p3lo[iso], v1 = p3hi[iso];
#ifdef SKIP_MUL
            add192(acc, f + x0 + x1 + v0 + v1, 0, 0);   // timing experiment: no 192-bit product
#else
            u64 p0 = x0 * v0, h0 = __umul64hi(x0, v0);
            u64 m2 = x1 * v0, m2h = __umul64hi(x1, v0);
            if (v1 == 0) {                  // 3^iso < 2^64 (iso <= 40): the usual case
                u64 s1 = h0 + m2;
                add192(acc, p0, s1, m2h + (s1 < m2));
            } else {
                u64 m1 = x0 * v1, m1h = __umul64hi(x0, v1);
                u64 s1 = h0 + m1; u64 c1 = s1 < m1;
                u64 s2b = s1 + m2; c1 += s2b < m2;
                add192(acc, p0, s2b, m1h + m2h + x1 * v1 + c1);
            }
#endif
        }
        __syncwarp();
    }

    // block reduction: warp shuffles, then one partial per warp
    for (int off = 16; off; off >>= 1) {
        u64 q0 = __shfl_down_sync(0xFFFFFFFF, acc.w[0], off);
        u64 q1 = __shfl_down_sync(0xFFFFFFFF, acc.w[1], off);
        u64 q2 = __shfl_down_sync(0xFFFFFFFF, acc.w[2], off);
        add192(acc, q0, q1, q2);
    }
    if (lane == 0) { red[0][warp] = acc.w[0]; red[1][warp] = acc.w[1]; red[2][warp] = acc.w[2]; }
    __syncthreads();
    if (tid == 0) {
        U192 t = {{0, 0, 0}};
        for (int i = 0; i < NWARP; i++) add192(t, red[0][i], red[1][i], red[2][i]);
        // times 3^base_iso (< 2^26), then << base_c2 (<= 32); the result is < 2^133
        const u64 W = p3lo[s_base_iso];
        u64 r0 = t.w[0] * W, c = __umul64hi(t.w[0], W);
        u64 r1 = t.w[1] * W; u64 c1 = __umul64hi(t.w[1], W) + ((r1 += c) < c);
        u64 r2 = t.w[2] * W + c1;
        const u32 sh = s_base_c2;
        if (sh) { r2 = (r2 << sh) | (r1 >> (64 - sh)); r1 = (r1 << sh) | (r0 >> (64 - sh)); r0 <<= sh; }
        t.w[0] = r0; t.w[1] = r1; t.w[2] = r2;
        u64 *o = out + 3 * ((u64)blockIdx.y * gridDim.x + blockIdx.x);
        o[0] = t.w[0]; o[1] = t.w[1]; o[2] = t.w[2];
    }
}

// ------------------------------------------------------------------------------------------------
// Host
// ------------------------------------------------------------------------------------------------
static std::string dec192(const U192 &v) {
    u64 w[3] = {v.w[0], v.w[1], v.w[2]};
    std::vector<u64> parts;
    const u64 B = 10000000000000000000ULL;
    while (w[0] || w[1] || w[2]) {
        unsigned __int128 rem = 0;
        for (int i = 2; i >= 0; i--) {
            unsigned __int128 cur = (rem << 64) | w[i];
            w[i] = (u64)(cur / B);
            rem = cur % B;
        }
        parts.push_back((u64)rem);
    }
    if (parts.empty()) return "0";
    std::string s = std::to_string(parts.back());
    char buf[32];
    for (int i = (int)parts.size() - 2; i >= 0; i--) { snprintf(buf, sizeof buf, "%019llu", (unsigned long long)parts[i]); s += buf; }
    return s;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: a8 REPS.bin --from I --to J [--batch B] [--out FILE] | --check REF.txt\n"); return 1; }
    auto opt = [&](const char *name, const char *def) -> std::string {
        for (int i = 2; i + 1 < argc; i++) if (!strcmp(argv[i], name)) return argv[i + 1];
        return def ? def : "";
    };
    // reps file: (u32 rep, u32 orbit size) pairs
    std::vector<u32> rep, wt;
    {
        FILE *fp = fopen(argv[1], "rb");
        if (!fp) { perror(argv[1]); return 1; }
        u32 pr[2];
        while (fread(pr, 4, 2, fp) == 2) { rep.push_back(pr[0]); wt.push_back(pr[1]); }
        fclose(fp);
    }
    fprintf(stderr, "%zu reps loaded\n", rep.size());

    // tables for Q_4
    u16 n4[16];
    for (int p = 0; p < 16; p++) { n4[p] = 0; for (int j = 0; j < 4; j++) n4[p] |= 1 << (p ^ (1 << j)); }
    auto N4 = [&](u32 y) { u32 r = 0; for (int p = 0; p < 16; p++) if (y >> p & 1) r |= n4[p]; return r; };
    std::vector<u16> i4;
    for (u32 s = 0; s < 65536; s++) if ((s & N4(s)) == 0) i4.push_back(s);
    if ((int)i4.size() != NI4) { fprintf(stderr, "bad I_4 count %zu\n", i4.size()); return 1; }
    std::vector<u16> fm(65536, 0);
    { std::vector<u32> f(65536, 0); for (u16 s : i4) f[s] = 1;
      for (int b = 0; b < 16; b++) for (u32 x = 0; x < 65536; x++) if (x >> b & 1) f[x] += f[x ^ (1u << b)];
      for (u32 x = 0; x < 65536; x++) fm[x] = (u16)f[x]; }
    std::vector<u16> i3;   // independent sets of Q_3 = the bytes that occur in I_4 sets
    for (u32 s = 0; s < 256; s++) {   // Q_3 on positions 0..7
        u32 n = 0; for (int p = 0; p < 8; p++) if (s >> p & 1) for (int j = 0; j < 3; j++) n |= 1u << (p ^ (1 << j));
        if ((s & n) == 0) i3.push_back(s);
    }
    if ((int)i3.size() != NI3) { fprintf(stderr, "bad I_3 count %zu\n", i3.size()); return 1; }
    std::vector<u8> hid(NI4), lid(NI4);
    for (int i = 0; i < NI4; i++) {
        int h = std::find(i3.begin(), i3.end(), i4[i] >> 8) - i3.begin();
        if (h == NI3 || std::find(i3.begin(), i3.end(), i4[i] & 255) == i3.end()) { fprintf(stderr, "I_4 set not a pair of I_3 sets\n"); return 1; }
        hid[i] = h;
        lid[i] = std::find(i3.begin(), i3.end(), i4[i] & 255) - i3.begin();
    }
    u64 hmask[256];
    for (u32 H = 0; H < 256; H++) { hmask[H] = 0; for (int h = 0; h < NI3; h++) if ((i3[h] & ~H) == 0) hmask[H] |= 1ULL << h; }
    std::vector<u32> lstoff(65537, 0); std::vector<u16> lst;
    for (u32 x = 0; x < 65536; x++) { for (u16 s : i4) if ((s & ~x) == 0) lst.push_back(s); lstoff[x + 1] = lst.size(); }
    u64 p3lo[65], p3hi[65];
    { unsigned __int128 v = 1; for (int i = 0; i <= 64; i++) { p3lo[i] = (u64)v; p3hi[i] = (u64)(v >> 64); v *= 3; } }

    CK(cudaMemcpyToSymbol(c_p3lo, p3lo, sizeof p3lo));
    CK(cudaMemcpyToSymbol(c_p3hi, p3hi, sizeof p3hi));
    CK(cudaMemcpyToSymbol(c_n4, n4, sizeof n4));
    CK(cudaMemcpyToSymbol(c_hmask, hmask, sizeof hmask));
    { u8 bo[256]; for (int i = 0; i < 256; i++) bo[i] = i;
      std::stable_sort(bo, bo + 256, [](u8 x, u8 y) { return __builtin_popcount(x) < __builtin_popcount(y); });
      CK(cudaMemcpyToSymbol(c_border, bo, sizeof bo)); }
    { u8 b[NI3]; for (int i = 0; i < NI3; i++) b[i] = (u8)i3[i]; CK(cudaMemcpyToSymbol(c_i3, b, sizeof b));
      u8 ix[256]; memset(ix, 255, sizeof ix); for (int i = 0; i < NI3; i++) ix[i3[i]] = i; CK(cudaMemcpyToSymbol(c_i3idx, ix, sizeof ix)); }
    { u8 ao[256]; for (int i = 0; i < 256; i++) ao[i] = i;
      std::stable_sort(ao, ao + 256, [](u8 x, u8 y) { return __builtin_popcount(x) > __builtin_popcount(y); });
      CK(cudaMemcpyToSymbol(c_aorder, ao, sizeof ao)); }
    u32 *d_reps; u16 *d_fm, *d_im, *d_lst; u32 *d_lstoff; u64 *d_out; u8 *d_hid;
    CK(cudaMalloc(&d_hid, NI4)); CK(cudaMemcpy(d_hid, hid.data(), NI4, cudaMemcpyHostToDevice));
    u8 *d_lid; CK(cudaMalloc(&d_lid, NI4)); CK(cudaMemcpy(d_lid, lid.data(), NI4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_reps, rep.size() * 4)); CK(cudaMemcpy(d_reps, rep.data(), rep.size() * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_fm, 65536 * 2)); CK(cudaMemcpy(d_fm, fm.data(), 65536 * 2, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_im, NI4 * 2)); CK(cudaMemcpy(d_im, i4.data(), NI4 * 2, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_lstoff, 65537 * 4)); CK(cudaMemcpy(d_lstoff, lstoff.data(), 65537 * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_lst, lst.size() * 2)); CK(cudaMemcpy(d_lst, lst.data(), lst.size() * 2, cudaMemcpyHostToDevice));

    // which rep indices to do
    std::vector<int> todo;
    std::vector<std::string> expect;
    std::string check = opt("--check", nullptr), outpath = opt("--out", nullptr);
    if (!check.empty()) {
        FILE *fp = fopen(check.c_str(), "r");
        if (!fp) { perror(check.c_str()); return 1; }
        char line[512];
        while (fgets(line, sizeof line, fp)) {
            int i; unsigned rr, ww; char s[256];
            if (sscanf(line, "%d %u %u %255s", &i, &rr, &ww, s) == 4) { todo.push_back(i); expect.push_back(s); }
        }
        fclose(fp);
    } else {
        int from = atoi(opt("--from", "0").c_str()), to = atoi(opt("--to", std::to_string(rep.size()).c_str()).c_str());
        std::set<int> done;
        if (!outpath.empty()) {
            FILE *fp = fopen(outpath.c_str(), "r");
            // drop a partial last line (cut off by a kill) so that nothing is ever appended to it
            if (fp) {
                std::string all; char buf[1 << 16]; size_t r;
                while ((r = fread(buf, 1, sizeof buf, fp)) > 0) all.append(buf, r);
                size_t nl = all.rfind('\n');
                size_t keep = nl == std::string::npos ? 0 : nl + 1;
                if (keep != all.size()) {
                    fprintf(stderr, "dropping a partial last line (%zu bytes) from %s\n", all.size() - keep, outpath.c_str());
                    if (truncate(outpath.c_str(), keep) != 0) { perror("truncate"); return 1; }
                }
                rewind(fp);
            }
            // only complete lines count
            if (fp) { char line[512]; while (fgets(line, sizeof line, fp)) {
                int i; unsigned rr, ww; char sum[256];
                size_t L = strlen(line);
                if (L && line[L - 1] == '\n' && sscanf(line, "%d %u %u %255s", &i, &rr, &ww, sum) == 4) done.insert(i);
            } fclose(fp); }
        }
        int stride = atoi(opt("--stride", "1").c_str());   // e.g. for timing a spread-out sample
        for (int i = from; i < to && i < (int)rep.size(); i += stride) if (!done.count(i)) todo.push_back(i);
        fprintf(stderr, "%zu reps to do (%zu already in %s)\n", todo.size(), done.size(), outpath.c_str());
    }
    int batch = atoi(opt("--batch", "2").c_str());
    int sym = 0;
    for (int i = 2; i < argc; i++) if (!strcmp(argv[i], "--sym")) sym = 1;
    fprintf(stderr, "swap symmetry: %s\n", sym ? "on (weighted rep sums)" : "off (plain rep sums)");
    FILE *fo = outpath.empty() ? stdout : fopen(outpath.c_str(), "a");
    CK(cudaMalloc(&d_out, (size_t)batch * 65536 * 3 * 8));
    std::vector<u64> h_out((size_t)batch * 65536 * 3);
    // the kernel takes a contiguous range of rep indices, so copy the batch's reps into a small buffer
    u32 *d_brep; CK(cudaMalloc(&d_brep, batch * 4));
    int bad = 0;
    // stall watchdog: a batch normally takes well under a second; if the GPU has not finished it
    // within STALL seconds, give up (exit code 3) so that the chunk loop restarts us (resuming)
    const double STALL = atof(opt("--stall", "120").c_str());
    cudaEvent_t ev; CK(cudaEventCreateWithFlags(&ev, cudaEventDisableTiming));
    auto t0 = std::chrono::steady_clock::now();
    for (size_t k = 0; k < todo.size(); k += batch) {
        int nb = std::min((size_t)batch, todo.size() - k);
        std::vector<u32> br(nb);
        for (int j = 0; j < nb; j++) br[j] = rep[todo[k + j]];
        CK(cudaMemcpy(d_brep, br.data(), nb * 4, cudaMemcpyHostToDevice));
        dim3 grid(65536, nb);
        slice_kernel<<<grid, THREADS>>>(d_brep, 0, d_fm, d_im, d_hid, d_lid, d_lstoff, d_lst, d_out, sym);
        CK(cudaGetLastError());
        CK(cudaEventRecord(ev));
        auto tl = std::chrono::steady_clock::now();
        while (true) {
            cudaError_t q = cudaEventQuery(ev);
            if (q == cudaSuccess) break;
            if (q != cudaErrorNotReady) { fprintf(stderr, "CUDA error %s while waiting\n", cudaGetErrorString(q)); fflush(stderr); _exit(2); }
            if (std::chrono::duration<double>(std::chrono::steady_clock::now() - tl).count() > STALL) {
                fprintf(stderr, "STALL: batch starting at rep %d not finished after %.0f s; exiting for a restart\n", todo[k], STALL);
                fflush(stderr);
                _exit(3);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
        }
        CK(cudaMemcpy(h_out.data(), d_out, (size_t)nb * 65536 * 3 * 8, cudaMemcpyDeviceToHost));
        std::string lines;   // the batch's result lines, written in one go
        for (int j = 0; j < nb; j++) {
            U192 s = {{0, 0, 0}};
            for (int y3 = 0; y3 < 65536; y3++) { u64 *o = &h_out[3 * ((size_t)j * 65536 + y3)]; add192(s, o[0], o[1], o[2]); }
            int i = todo[k + j];
            std::string d = dec192(s);
            if (!check.empty()) {
                bool ok = d == expect[k + j];
                bad += !ok;
                fprintf(stderr, "rep #%d: gpu %s  cpu %s  %s\n", i, d.c_str(), expect[k + j].c_str(), ok ? "OK" : "MISMATCH");
            } else {
                char buf[64]; snprintf(buf, sizeof buf, "%d %u %u ", i, rep[i], wt[i]);
                lines += buf + d + "\n";
            }
        }
        if (!lines.empty()) { fwrite(lines.data(), 1, lines.size(), fo); fflush(fo); }
        double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        size_t done = k + nb;
        if (check.empty() && (done % (batch * 50) == 0 || done == todo.size()))
            fprintf(stderr, "[%8.0fs] %zu/%zu reps  (%.3f s/rep, ETA %.1f h)\n", el, done, todo.size(), el / done, el / done * (todo.size() - done) / 3600);
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "done: %zu reps in %.1fs (%.3f s/rep)%s\n", todo.size(), el, el / std::max<size_t>(1, todo.size()),
            check.empty() ? "" : (bad ? "  CHECK FAILED" : "  ALL MATCH"));
    return bad ? 1 : 0;
}
