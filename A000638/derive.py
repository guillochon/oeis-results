"""Derive A000638(21) from A000637(21).

    python derive.py

A000638(n) = Sum_{k=0..n} A000637(k): a permutation group of degree n is a fixed-point-free group on its set of
moved points (see README.md). The relation is checked on every known term before the new one is printed.
"""

# Fixed-point-free permutation groups of degree n, n = 0..21 (OEIS DATA, Oct 2026). a(21) = 140399143 is
# Jeffrey Ketchersid's S21 computation (added Oct 12 2026).
A000637 = [1, 0, 1, 2, 7, 8, 37, 40, 200, 258, 1039, 1501, 7629, 10109, 54322, 83975, 527036, 780193, 5808293,
           9470582, 88249363, 140399143]

# Permutation groups of degree n, n = 0..20 (OEIS DATA, Oct 2026), used as a check.
A000638_KNOWN = [1, 1, 2, 4, 11, 19, 56, 96, 296, 554, 1593, 3094, 10723, 20832, 75154, 159129, 686165, 1466358,
                 7274651, 16745233, 104994596]


def main():
    partial, s = [], 0
    for x in A000637:
        s += x
        partial.append(s)
    n = len(A000638_KNOWN)
    assert partial[:n] == A000638_KNOWN
    print(f"partial sums of A000637 reproduce all {n} known terms A000638(0..{n - 1})")
    print(f"A000638({n}) = {A000638_KNOWN[-1]} + {A000637[n]} = {partial[n]}")


if __name__ == "__main__":
    main()
