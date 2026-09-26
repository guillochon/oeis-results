"""Write the b-file b236554.txt (n = 1..500) from the proved formula a(n) = 2^(2n+2) + 32 (n >= 3),
after checking it against the direct count of count.py for n <= 13.
Usage: python bfile.py
"""
from pathlib import Path

from count import count

HERE = Path(__file__).resolve().parent


def a(n):
    return {1: 8, 2: 64}.get(n, 2 ** (2 * n + 2) + 32)


if __name__ == "__main__":
    assert all(a(n) == count(n) for n in range(1, 14)), "formula disagrees with direct count"
    (HERE / "b236554.txt").write_text("".join(f"{n} {a(n)}\n" for n in range(1, 501)))
    print("wrote b236554.txt (n = 1..500); formula checked against direct count for n <= 13")
