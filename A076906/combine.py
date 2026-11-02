"""Combine the transfer-matrix runs into the six sequences of the family.

Inputs (results/, written by run_all.sh and the fix / first runs in README):
  q1_N, q4_N+1, q3_N          paths from Q3 ending in Q1 / Q4 / Q3
  q2q4_N, q2q1_N+1, q2q2_N    reflected cross-checks (start in Q2)
  fixax_q1_20, fixax_q2_20    half curves for axis-symmetric curves
  fixdg_q3_20, fixdg_q2_20, fixdg_q1_20   half curves for diagonal-symmetric curves
  q3first_M                   Q3 -> Q3 with first-crossing flags (A077556, even n)
Output: results/counts.json and b-files in bfiles/.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"


def load(pattern):
    """{n: [count, xfirst, xx?]} from the newest file matching pattern."""
    fs = sorted(R.glob(pattern), key=lambda p: int(re.findall(r"_(\d+)\.txt$", p.name)[0]))
    if not fs:
        return {}
    out = {}
    for line in fs[-1].read_text().splitlines():
        f = line.split()
        nums = []
        for x in f[1:]:
            if not x.isdigit():
                break
            nums.append(int(x))
        out[int(f[0])] = nums
    return out


def col(d, j=0):
    return {n: v[j] for n, v in d.items()}


def main():
    q1, q4, q3 = col(load("q1_*.txt")), col(load("q4_*.txt")), col(load("q3_*.txt"))
    x14, x41, x22 = col(load("q2q4_*.txt")), col(load("q2q1_*.txt")), col(load("q2q2_*.txt"))
    # reflection in the x axis: Q3->Q1 ~ Q2->Q4, Q3->Q4 ~ Q2->Q1, Q3->Q3 ~ Q2->Q2
    for a, b, name in ((q1, x14, "q1/q2q4"), (q4, x41, "q4/q2q1"), (q3, x22, "q3/q2q2")):
        common = sorted(set(a) & set(b))
        bad = [n for n in common if a[n] != b[n]]
        assert not bad, (name, bad)
        print(f"cross-check {name}: n <= {max(common)} agree")
    N = 0                                 # even n need the Q1 and Q3 runs, odd n the Q4 run
    while (N + 1 in q4) if (N + 1) % 2 else (N + 1 in q1 and N + 1 in q3):
        N += 1
    fa1, fa2 = col(load("fixax_q1_*.txt")), col(load("fixax_q2_*.txt"))
    fd3, fd2, fd1 = col(load("fixdg_q3_*.txt")), col(load("fixdg_q2_*.txt")), col(load("fixdg_q1_*.txt"))
    A005316 = {k: fa1[k] + fa2[k] for k in fa1}          # identity checked below
    A076876 = {k: fd3[k] + fd2[k] + fd1[k] for k in fd1}
    A005316[0] = A076876[0] = 1      # zero crossings: the engine does not count the empty curve
    # the half-curve counts are known sequences (b-files of A005316 and A076876)
    known5316 = [1, 1, 1, 2, 3, 8, 14, 42, 81, 262, 538, 1828, 3926, 13820, 30694, 110954, 252939,
                 933458, 2172830, 8152860, 19304190]
    known6876 = [1, 1, 2, 3, 8, 14, 43, 81, 272, 538, 1920, 3926, 14649, 30694, 118489, 252939, 1002994,
                 2172830, 8805410, 19304190, 79648888]
    assert [A005316[k] for k in range(21)] == known5316
    assert [A076876[k] for k in range(21)] == known6876
    first = load("q3first_*.txt")

    seqs = {k: {} for k in ("A076906", "A076907", "A077551", "A076875", "A077550", "A077556")}
    seqs["A076906"][0], seqs["A076907"][0], seqs["A077551"][0] = 0, 2, 0
    seqs["A076875"][0], seqs["A077550"][0], seqs["A077556"][0] = 1, 1, 0
    for n in range(1, N + 1):
        if n % 2:
            a906, a907, a875 = q4[n], 2 * q4[n], 2 * q4[n]
            a551 = a906                         # reflection in y = x
            a550 = (a907 + 2 * A005316[(n - 1) // 2]) // 4
            assert (a907 + 2 * A005316[(n - 1) // 2]) % 4 == 0
            a556 = a906
        else:
            a906, a907 = q1[n], q1[n] + q3[n]
            a875 = q1[n] + q3[n] // 2
            a551 = a907 // 2
            a550 = (a907 + 2 * A076876[n // 2]) // 4
            assert (a907 + 2 * A076876[n // 2]) % 4 == 0
            a556 = q1[n] // 2 + (q3[n] - first[n][2]) // 2 if n in first else None
        for k, v in zip(seqs, (a906, a907, a551, a875, a550, a556)):
            if v is not None:
                seqs[k][n] = v
    for k, v in seqs.items():                # keep terms contiguous from n = 0
        m = 0
        while m + 1 in v:
            m += 1
        seqs[k] = {n: v[n] for n in range(m + 1)}
    out = {k: [v[n] for n in sorted(v)] for k, v in seqs.items()}
    out["A005316_from_half_curves"] = [A005316[k] for k in sorted(A005316)]
    out["A076876_from_half_curves"] = [A076876[k] for k in sorted(A076876)]
    (R / "counts.json").write_text(json.dumps(out, indent=1))
    bdir = HERE / "bfiles"
    bdir.mkdir(exist_ok=True)
    for k in seqs:
        lines = [f"{n} {seqs[k][n]}" for n in sorted(seqs[k])]
        (bdir / f"b{k[1:]}.txt").write_text("\n".join(lines) + "\n")
    for k, v in out.items():
        print(k, len(v), v[:12], "...", v[-1])


if __name__ == "__main__":
    main()
