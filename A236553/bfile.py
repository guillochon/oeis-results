"""Write the b-file b236553.txt (n = 1..10000) from the proved closed form, after checking the
closed form against a direct count for every n <= 256 (verify.py).
Usage: python bfile.py
"""
from pathlib import Path

from sympy import factorint

from verify import count, local_A236553

HERE = Path(__file__).resolve().parent


def a(n):
    out = 1
    for p, e in factorint(n).items():
        out *= local_A236553(p, e)
    return out


if __name__ == "__main__":
    assert all(a(n) == count(n, (-1, 1, 1)) for n in range(1, 257)), "closed form disagrees"
    (HERE / "b236553.txt").write_text("".join(f"{n} {a(n)}\n" for n in range(1, 10001)))
    print("wrote b236553.txt (n = 1..10000); closed form checked against direct count for n <= 256")
