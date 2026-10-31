"""Derive A005226(19..21) as the inverse Euler transform of A000638(0..21).

    python derive.py

A000638 (molecular species of degree n) is the Euler transform of A005226 (atomic species of degree n), because
every molecular species is uniquely a product of atomic ones (see README.md). The inversion is checked on every
known term, and the Euler transform of the result is checked to give back A000638, before the new terms are
printed.
"""

# Fixed-point-free permutation groups of degree n, n = 0..21 (OEIS DATA, Oct 2026). a(19), a(20), a(21) are
# Jeffrey Ketchersid's S19, S20, S21 computations.
A000637 = [1, 0, 1, 2, 7, 8, 37, 40, 200, 258, 1039, 1501, 7629, 10109, 54322, 83975, 527036, 780193, 5808293,
           9470582, 88249363, 140399143]

# Permutation groups (= molecular species) of degree n, n = 0..20 (OEIS DATA, Oct 2026).
A000638_KNOWN = [1, 1, 2, 4, 11, 19, 56, 96, 296, 554, 1593, 3094, 10723, 20832, 75154, 159129, 686165, 1466358,
                 7274651, 16745233, 104994596]

# Atomic species of degree n, n = 0..18 (OEIS DATA, Oct 2026), used as a check.
A005226_KNOWN = [0, 1, 1, 2, 6, 6, 27, 20, 130, 124, 598, 641, 4850, 4772, 35625, 46074, 389839, 487408, 4617554]


def euler(b, N):
    """a with 1 + Sum_{n>=1} a(n) x^n = Prod_{k>=1} (1 - x^k)^(-b(k)), for n = 0..N."""
    c = [0] + [sum(d * b[d] for d in range(1, n + 1) if n % d == 0) for n in range(1, N + 1)]
    a = [1] + [0] * N
    for n in range(1, N + 1):
        s = sum(c[k] * a[n - k] for k in range(1, n + 1))
        assert s % n == 0
        a[n] = s // n
    return a


def inverse_euler(a):
    """b with euler(b) = a (a(0) = 1), b(0) = 0."""
    N = len(a) - 1
    c = [0] * (N + 1)
    for n in range(1, N + 1):
        c[n] = n * a[n] - sum(c[k] * a[n - k] for k in range(1, n))
    b = [0] * (N + 1)
    for n in range(1, N + 1):
        s = sum(d * b[d] for d in range(1, n) if n % d == 0)
        assert (c[n] - s) % n == 0
        b[n] = (c[n] - s) // n
    return b


def main():
    # A000638(21) = A000638(20) + A000637(21): partial sums of A000637 (a group of degree n is a fixed-point-free
    # group on its moved points).
    partial = [sum(A000637[:n + 1]) for n in range(len(A000637))]
    assert partial[:len(A000638_KNOWN)] == A000638_KNOWN
    A000638 = partial
    print(f"A000638(0..{len(A000638) - 1}) = partial sums of A000637; A000638(21) = {A000638[21]}")

    b = inverse_euler(A000638)
    k = len(A005226_KNOWN)
    assert b[:k] == A005226_KNOWN
    print(f"inverse Euler transform reproduces all {k} known terms A005226(0..{k - 1})")
    assert euler(b, len(b) - 1) == A000638
    print(f"Euler transform of A005226(0..{len(b) - 1}) gives back A000638(0..{len(b) - 1})")
    for n in range(k, len(b)):
        print(f"A005226({n}) = {b[n]}")


if __name__ == "__main__":
    main()
