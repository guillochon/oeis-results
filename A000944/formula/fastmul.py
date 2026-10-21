"""Exact truncated multiplication of multivariate series mod p (p < 2^26) by floating-point FFT.

mulmod(a, b, p, shape) = (a * b mod p) truncated to `shape`, for int64 arrays with entries in [0, p).
`shape` covers the trailing (series) axes; any leading axes are batch axes, multiplied pointwise with
numpy broadcasting (an operand with fewer dimensions gets leading axes of size 1). Each operand is split
into k limbs of `bits` bits (a = sum a_i 2^(bits i)); the limb products are exact in float64 as long as
the convolution sums stay well below 2^53, which `_split` guarantees from the number of terms, and every
result is checked to be within 0.25 of an integer. Operands are first cut to the bounding box of their
nonzero entries, and sparse operands (monomials, a few terms, possibly varying over the batch) are
handled by shifted adds.
"""
import numpy as np

try:
    from scipy import fft as _fft
    _KW = {'workers': 1}
except ImportError:            # numpy's pocketfft: same interface, no workers argument
    _fft = np.fft
    _KW = {}

SPARSE = 8                     # operands with at most this many nonzero terms are multiplied directly


def _fast_len(n):
    if hasattr(_fft, 'next_fast_len'):
        return _fft.next_fast_len(n, real=True)
    return 1 << (n - 1).bit_length()


def _bbox(a, nb):
    """(low, high) corners of the nonzero entries over the series axes, or None if a = 0."""
    m = a.any(axis=tuple(range(nb))) if nb else a
    nz = np.nonzero(m)
    if not len(nz[0]):
        return None
    return tuple(int(ix.min()) for ix in nz), tuple(int(ix.max()) + 1 for ix in nz)


def _split(nterms, pbits):
    """(limbs, bits): limb products summed over nterms stay below 2^(44)."""
    for k in (1, 2, 3, 4):
        bits = -(-pbits // k)
        if 2 * bits + nterms.bit_length() <= 44:
            return k, bits
    raise ValueError('series too large for float64 FFT multiplication')


def _sparse_mul(a, b, p, shape, bshape):
    """b has few nonzero positions (over all batch entries): sum of shifted copies of a."""
    nb = len(bshape)
    out = np.zeros(bshape + shape, dtype=np.int64)
    m = b.any(axis=tuple(range(nb))) if nb else b
    for idx in zip(*np.nonzero(m)):
        c = b[(Ellipsis,) + idx].reshape(b.shape[:nb] + (1,) * len(shape))
        n = [min(k, s - i) for k, i, s in zip(a.shape[nb:], idx, shape)]
        if min(n) <= 0:
            continue
        dst = (Ellipsis,) + tuple(slice(i, i + k) for i, k in zip(idx, n))
        src = (Ellipsis,) + tuple(slice(0, k) for k in n)
        out[dst] = (out[dst] + a[src] * c) % p
    return out


def mulmod(a, b, p, shape=None):
    shape = tuple(shape or a.shape[-max(a.ndim, b.ndim):])
    d = len(shape)
    nd = max(a.ndim, b.ndim)
    a = a.reshape((1,) * (nd - a.ndim) + a.shape)
    b = b.reshape((1,) * (nd - b.ndim) + b.shape)
    nb = nd - d
    bshape = tuple(np.broadcast_shapes(a.shape[:nb], b.shape[:nb]))
    ba, bb = _bbox(a, nb), _bbox(b, nb)
    if ba is None or bb is None:
        return np.zeros(bshape + shape, dtype=np.int64)
    # cut operands to their bounding boxes (low and high ends), and to the output shape
    lo = tuple(x + y for x, y in zip(ba[0], bb[0]))       # lowest corner of the product
    if any(l >= s for l, s in zip(lo, shape)):
        return np.zeros(bshape + shape, dtype=np.int64)
    a = a[(Ellipsis,) + tuple(slice(l, min(h, s - m)) for l, h, s, m in zip(*ba, shape, bb[0]))]
    b = b[(Ellipsis,) + tuple(slice(l, min(h, s - m)) for l, h, s, m in zip(*bb, shape, ba[0]))]
    sub = tuple(s - l for s, l in zip(shape, lo))           # product window, relative to lo
    res = np.zeros(bshape + shape, dtype=np.int64)
    res[(Ellipsis,) + tuple(slice(l, None) for l in lo)] = _mul_window(a, b, p, sub, nb, bshape)
    return res


def _mul_window(a, b, p, shape, nb, bshape):
    """a * b truncated to shape, for operands already cut to their bounding boxes."""
    ba, bb = a.shape[nb:], b.shape[nb:]
    if a.min() < 0 or a.max() >= p:      # callers may pass unreduced sums or negatives
        a = a % p
    if b.min() < 0 or b.max() >= p:
        b = b % p
    if np.count_nonzero(b.any(axis=tuple(range(nb))) if nb else b) <= SPARSE:
        return _sparse_mul(a, b, p, shape, bshape)
    if np.count_nonzero(a.any(axis=tuple(range(nb))) if nb else a) <= SPARSE:
        return _sparse_mul(b, a, p, shape, bshape)
    full = tuple(min(x + y - 1, s) for x, y, s in zip(ba, bb, shape))   # the part of the product we keep
    fft_shape = tuple(_fast_len(x + y - 1) for x, y in zip(ba, bb))
    nterms = int(np.prod([min(x, y) for x, y in zip(ba, bb)]))
    k, bits = _split(nterms, max(1, (p - 1).bit_length()))
    mask = (1 << bits) - 1
    axes = tuple(range(nb, a.ndim))
    FA = [_fft.rfftn(((a >> (bits * i)) & mask).astype(np.float64), fft_shape, axes=axes, **_KW) for i in range(k)]
    FB = [_fft.rfftn(((b >> (bits * i)) & mask).astype(np.float64), fft_shape, axes=axes, **_KW) for i in range(k)]
    keep = (Ellipsis,) + tuple(slice(0, x) for x in full)
    out = np.zeros(bshape + full, dtype=np.int64)
    for s in range(2 * k - 1):
        F = sum(FA[i] * FB[s - i] for i in range(max(0, s - k + 1), min(s, k - 1) + 1))
        c = _fft.irfftn(F, fft_shape, axes=axes, **_KW)[keep]
        r = np.rint(c)
        err = np.abs(c - r).max() if c.size else 0.0
        if err > 0.25:
            raise ArithmeticError(f'FFT rounding error {err:.3f}')
        out = (out + (r.astype(np.int64) % p) * pow(2, bits * s, p)) % p
    res = np.zeros(bshape + shape, dtype=np.int64)
    res[keep] = out
    return res
