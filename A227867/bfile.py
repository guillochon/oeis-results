"""Write the b-file b227867.txt (n = 1..10000) from the proved closed form, after checking the
closed form against a direct count for every n <= 256 (verify.py, shared with A236553).
Usage: python bfile.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "A236553")]

from sympy import factorint  # noqa: E402

from verify import count, local_A227867  # noqa: E402


def a(n):
    out = 1
    for p, e in factorint(n).items():
        out *= local_A227867(p, e)
    return out


if __name__ == "__main__":
    assert all(a(n) == count(n, (-1, -1, -1)) for n in range(1, 257)), "closed form disagrees"
    (HERE / "b227867.txt").write_text("".join(f"{n} {a(n)}\n" for n in range(1, 10001)))
    print("wrote b227867.txt (n = 1..10000); closed form checked against direct count for n <= 256")
