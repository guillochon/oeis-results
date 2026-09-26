"""Explicit construction of long A157615 paths by frame induction: optimal for odd n, and of the
conjectured length n^2 - n + 2 for even n.

ODD n. Invariant: P(n) is a valid alternating path on the n X n board of length n^2 - 2n + 4 whose last
two squares are (n-1, 0) -> (n-1, 1) (it ends on the right edge, second row, arriving downward).

Step n -> n+2: append the "frame route" through the new columns n, n+1 and rows n, n+1:
    (n, 1);
    rows r = 2..n of the right strip: (n, r), (n+1, r) if r is even, else (n+1, r), (n, r);
    (n, n+1);
    then leftwards along the bottom strip, column c = n-1..0:
        (c, n+1), (c, n) if n-1-c is even, else (c, n), (c, n+1).
This adds exactly 4n squares and ends (0, n+1) -> (0, n). Rotating the (n+2)-board by 180 degrees
turns that ending into (n+1, 0) -> (n+1, 1), so the invariant holds for n + 2.

Together with the four-colouring upper bound, this proves a(n) = n^2 - 2n + 4 for all odd n >= 3.

EVEN n. Invariant: E(n) has length n^2 - n + 2 and ends (n-1, 1) -> (n-1, 0) (top-right corner,
arriving upward). Step n -> n+2 appends
    (n, 0);
    rows r = 1..n of the right strip: (n, r), (n+1, r) if r is odd, else (n+1, r), (n, r);
    (n, n+1);
    then leftwards along the bottom strip, column c = n-1..0:
        (c, n+1), (c, n) if n-1-c is even, else (c, n), (c, n+1).
This adds exactly 4n + 2 squares (only (n+1, 0) and (n+1, n+1) are missed) and ends
(0, n) -> (0, n+1); rotating by 180 degrees restores the invariant. So a(n) >= n^2 - n + 2 for all
even n >= 4, the conjectured value (the matching upper bound for even n is still open).

    python construct.py 301     # build and validate P(n), E(n) for n = 4..301
"""
import sys

from path import check, conjectured, draw

# Base cases, found by exhaustive search (brute.py); both end (n-1, 0) -> (n-1, 1).
BASE = {
    5: [(1, 0), (0, 0), (0, 1), (1, 1), (1, 2), (0, 2), (0, 3), (1, 3), (1, 4), (2, 4), (2, 3),
        (3, 3), (3, 2), (2, 2), (2, 1), (3, 1), (3, 0), (4, 0), (4, 1)],
}


BASE_EVEN = {
    4: [(3, 3), (3, 2), (2, 2), (2, 3), (1, 3), (1, 2), (0, 2), (0, 1), (1, 1), (1, 0), (2, 0), (2, 1),
        (3, 1), (3, 0)],
}


def frame_route_even(n):
    route = [(n, 0)]
    for r in range(1, n + 1):
        route += [(n, r), (n + 1, r)] if r % 2 == 1 else [(n + 1, r), (n, r)]
    route.append((n, n + 1))
    for c in range(n - 1, -1, -1):
        route += [(c, n + 1), (c, n)] if (n - 1 - c) % 2 == 0 else [(c, n), (c, n + 1)]
    return route


def frame_route(n):
    route = [(n, 1)]
    for r in range(2, n + 1):
        route += [(n, r), (n + 1, r)] if r % 2 == 0 else [(n + 1, r), (n, r)]
    route.append((n, n + 1))
    for c in range(n - 1, -1, -1):
        route += [(c, n + 1), (c, n)] if (n - 1 - c) % 2 == 0 else [(c, n), (c, n + 1)]
    return route


def rotate(N, path):
    return [(N - 1 - x, N - 1 - y) for x, y in path]


def P(n):
    """Optimal path for odd n >= 5 (satisfying the invariant)."""
    assert n % 2 == 1 and n >= 5
    p = BASE[5]
    for m in range(5, n, 2):
        p = rotate(m + 2, p + frame_route(m))
    return p


def E(n):
    """Path of length n^2 - n + 2 for even n >= 4 (satisfying the even invariant)."""
    assert n % 2 == 0 and n >= 4
    p = BASE_EVEN[4]
    for m in range(4, n, 2):
        p = rotate(m + 2, p + frame_route_even(m))
    return p


def invariant_ok(n, p):
    end = [(n - 1, 0), (n - 1, 1)] if n % 2 else [(n - 1, 1), (n - 1, 0)]
    return check(n, p) is None and len(p) == conjectured(n) and p[-2:] == end


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 101
    assert invariant_ok(5, BASE[5]) and invariant_ok(4, BASE_EVEN[4]), "base case"
    bad = [n for n in range(5, top + 1, 2) if not invariant_ok(n, P(n))]
    print(f"odd n = 5..{top}: " + ("every P(n) is valid, has length n^2 - 2n + 4, and satisfies the invariant"
                                  if not bad else f"FAILED for n = {bad[:10]}"))
    bad = [n for n in range(4, top + 1, 2) if not invariant_ok(n, E(n))]
    print(f"even n = 4..{top}: " + ("every E(n) is valid, has length n^2 - n + 2, and satisfies the invariant"
                                   if not bad else f"FAILED for n = {bad[:10]}"))
    for n in (7, 8):
        q = P(n) if n % 2 else E(n)
        print(f"\n{'P' if n % 2 else 'E'}({n}), length {len(q)}:\n{draw(n, q)}")
