"""GPU (CuPy) replacement for refl.S2: truncated 2D series mod p with a leading batch axis.

Multiplication is a direct truncated 2D convolution in a CUDA kernel (one block per batch entry,
both operands in shared memory, exact 64-bit accumulation: products < 2^48, at most X*Y <= 2^11
terms, then one reduction mod p).  For the small boxes of the half-map grammar (about 21 x 37) this
is far faster than FFT on a consumer GPU, whose float64 throughput is low.  Operands broadcast over
the batch axis (an unbatched operand has batch size 1).
"""
import cupy as cp
import numpy as np

_SRC = r'''
extern "C" __global__ void conv2(const long long* A, const long long* B, long long* C,
                                 int Xa, int Ya, int Xb, int Yb, int X, int Y,
                                 long long sa, long long sb, unsigned long long p)
{
    extern __shared__ unsigned long long sh[];
    unsigned long long* a = sh;
    unsigned long long* b = sh + Xa * Ya;
    const long long bt = blockIdx.x;
    const long long* Ag = A + bt * sa;
    const long long* Bg = B + bt * sb;
    for (int t = threadIdx.x; t < Xa * Ya; t += blockDim.x) a[t] = (unsigned long long)Ag[t];
    for (int t = threadIdx.x; t < Xb * Yb; t += blockDim.x) b[t] = (unsigned long long)Bg[t];
    __syncthreads();
    for (int o = threadIdx.x; o < X * Y; o += blockDim.x) {
        int i = o / Y, j = o % Y;
        unsigned long long acc = 0;
        int k0 = i - Xb + 1; if (k0 < 0) k0 = 0;
        int k1 = i < Xa - 1 ? i : Xa - 1;
        int l0 = j - Yb + 1; if (l0 < 0) l0 = 0;
        int l1 = j < Ya - 1 ? j : Ya - 1;
        for (int k = k0; k <= k1; k++) {
            const unsigned long long* ar = a + k * Ya;
            const unsigned long long* br = b + (i - k) * Yb + j;
            for (int l = l0; l <= l1; l++) acc += ar[l] * br[-l];
        }
        C[bt * (long long)(X * Y) + o] = (long long)(acc % p);
    }
}
'''
_conv2 = cp.RawKernel(_SRC, 'conv2')


def _as3(a, X, Y):
    """cupy int64 array reduced mod p elsewhere; returns (batch, x, y) view clipped to (X, Y)."""
    if a.ndim == 2:
        a = a[None]
    elif a.ndim > 3:
        a = a.reshape((-1,) + a.shape[-2:])
    return a[:, :X, :Y]


class GS2:
    """Drop-in for refl.S2 on the GPU (arrays are cupy int64, entries in [0, p))."""

    def __init__(self, VX, FY, p):
        self.VX, self.FY, self.p = VX, FY, p

    def zero(self):
        return cp.zeros((self.VX + 1, self.FY + 1), dtype=cp.int64)

    def mono(self, i, j, c=1):
        a = self.zero()
        if i <= self.VX and j <= self.FY:
            a[i, j] = c % self.p
        return a

    def mul(self, a, b):
        p = self.p
        X, Y = self.VX + 1, self.FY + 1
        a = cp.asarray(a)
        b = cp.asarray(b)
        lead = np.broadcast_shapes(a.shape[:-2], b.shape[:-2])
        a3 = cp.ascontiguousarray(_as3(a, X, Y) % p)
        b3 = cp.ascontiguousarray(_as3(b, X, Y) % p)
        Ba, Bb = a3.shape[0], b3.shape[0]
        Bt = max(Ba, Bb)
        assert Ba in (1, Bt) and Bb in (1, Bt)
        out = cp.empty((Bt, X, Y), dtype=cp.int64)
        Xa, Ya = a3.shape[1:]
        Xb, Yb = b3.shape[1:]
        shmem = 8 * (Xa * Ya + Xb * Yb)
        threads = min(1024, max(32, ((X * Y + 31) // 32) * 32))
        _conv2((Bt,), (threads,), (a3, b3, out, np.int32(Xa), np.int32(Ya), np.int32(Xb), np.int32(Yb),
                                    np.int32(X), np.int32(Y),
                                    np.int64(0 if Ba == 1 else Xa * Ya), np.int64(0 if Bb == 1 else Xb * Yb),
                                    np.uint64(p)), shared_mem=shmem)
        return out.reshape(tuple(lead) + (X, Y)) if lead else out.reshape(X, Y)


def selftest():
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'A000944', 'formula'))
    from refl import S2
    p = 16777213
    R, G = S2(20, 36, p), GS2(20, 36, p)
    rng = np.random.default_rng(1)
    for shape_a, shape_b in (((21, 37), (21, 37)), ((5, 21, 37), (21, 37)), ((5, 21, 37), (5, 21, 37)),
                             ((21, 37), (7, 21, 37))):
        a = rng.integers(0, p, shape_a)
        b = rng.integers(0, p, shape_b)
        want = R.mul(a, b)
        got = cp.asnumpy(G.mul(cp.asarray(a), cp.asarray(b)))
        assert (np.broadcast_to(want, got.shape) == got).all(), (shape_a, shape_b)
    print('gpu conv2 matches fastmul.mulmod')


if __name__ == '__main__':
    selftest()
