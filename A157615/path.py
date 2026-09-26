"""Validator and drawing helpers for A157615 paths.

A path is a list of (x, y) squares on the n X n board (0 <= x, y < n). It is valid iff all squares
are distinct and on the board, consecutive squares are unit steps apart, and the steps alternate
horizontal / vertical.
"""


def check(n, path):
    """Return None if `path` is a valid alternating self-avoiding path on n X n, else a reason."""
    if len(set(path)) != len(path):
        return "revisits a square"
    for x, y in path:
        if not (0 <= x < n and 0 <= y < n):
            return f"({x},{y}) is off the board"
    for i in range(len(path) - 1):
        (x1, y1), (x2, y2) = path[i], path[i + 1]
        dx, dy = abs(x2 - x1), abs(y2 - y1)
        if dx + dy != 1:
            return f"step {i} is not a unit step"
        if i > 0:
            px, _ = path[i - 1]
            prev_horizontal = px != x1
            if prev_horizontal == (dx == 1):
                return f"steps {i - 1} and {i} do not alternate"
    return None


def draw(n, path):
    g = [["." for _ in range(n)] for _ in range(n)]
    for i, (x, y) in enumerate(path):
        g[y][x] = "S" if i == 0 else "E" if i == len(path) - 1 else "#"
    return "\n".join(" ".join(r) for r in g)


def conjectured(n):
    """David Wilson's conjecture (the odd case is a proved upper bound, see README)."""
    if n == 1:
        return 1
    return n * n - n + 2 if n % 2 == 0 else n * n - 2 * n + 4
