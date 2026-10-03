"""Assemble b-files from the Rust residue outputs, checking them against the earlier Python values.

    python make_bfiles.py

  bfile_stacking.txt          stacking, distinct sizes, n = 0..   (from stack_residues.txt, 2 primes held back)
  bfile_repeat_nesting.txt    repeated sizes, nesting, n = 0..    (largest finished rnest_<N>.txt run)
  bfile_repeat_stacking.txt   repeated sizes, stacking, n = 0..   (largest finished rstack_<N>.txt run)
and refreshes repeat_nest_terms.txt / repeat_stack_terms.txt (n >= 1, one per line) used by the pages.
"""
import re
from pathlib import Path

from crt import reconstruct

HERE = Path(__file__).resolve().parent


def write_bfile(path, terms):
    path.write_text("".join(f"{n} {x}\n" for n, x in enumerate(terms)), encoding="utf-8")


def latest_run(prefix):
    """Largest N whose run printed every n = 1..N."""
    best = None
    for f in HERE.glob(f"{prefix}_*.txt"):
        m = re.fullmatch(rf"{prefix}_(\d+)\.txt", f.name)
        if not m:
            continue
        n = int(m.group(1))
        vals = reconstruct("repeat", f, 1)
        if sorted(vals) == list(range(1, n + 1)) and (best is None or n > best[0]):
            best = (n, vals)
    return best


def main():
    # stacking, distinct sizes
    old = [int(l.split()[1]) for l in (HERE / "bfile_stacking.txt").read_text().splitlines()]
    vals = reconstruct("stack", HERE / "stack_residues.txt", 2)
    terms = [1] + [vals[n] for n in range(1, max(vals) + 1)]
    assert terms[: len(old)] == old, "stacking: Rust disagrees with the Python b-file"
    write_bfile(HERE / "bfile_stacking.txt", terms)
    print(f"stacking: n = 0..{len(terms) - 1} (was 0..{len(old) - 1})")

    for model, prefix, termfile in (("nesting", "rnest", "repeat_nest_terms.txt"),
                                    ("stacking", "rstack", "repeat_stack_terms.txt")):
        n, vals = latest_run(prefix)
        terms = [1] + [vals[k] for k in range(1, n + 1)]
        py = [int(x) for x in (HERE / termfile).read_text().split()]
        assert terms[1: len(py) + 1] == py[: n], f"repeat {model}: Rust disagrees with Python"
        write_bfile(HERE / f"bfile_repeat_{model}.txt", terms)
        (HERE / termfile).write_text("".join(f"{x}\n" for x in terms[1:]), encoding="utf-8")
        print(f"repeated sizes, {model}: n = 0..{n} (Python had 1..{len(py)})")


if __name__ == "__main__":
    main()
