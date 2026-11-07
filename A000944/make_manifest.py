"""Write MANIFEST.md: inputs, residues and outputs of the a(50) run, with SHA-256 hashes.

    python attacks/A000944/make_manifest.py
"""
import hashlib
import platform
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
F = HERE / "formula"

DATA = ["b000944.txt", "b212438.txt", "b239893.txt", "b000287.txt", "b002840.txt", "b006445.txt"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def lines(p):
    return sum(1 for l in p.read_text().splitlines() if l.strip() and not l.startswith("#"))


def main():
    import numpy
    import scipy
    sys.path.insert(0, str(F))
    from ms import PRIMES
    from t3 import PRIMES24
    L = ["# A000944 run manifest", "",
         "Everything needed to replay the a(50) computation from this folder without network access, with",
         "SHA-256 hashes. Regenerate with `make_manifest.py` in the source repository.", "",
         "## Software", "",
         f"Python {platform.python_version()}, numpy {numpy.__version__}, scipy {scipy.__version__} "
         f"({platform.system()} {platform.release()}, i5-12600K). networkx and nauty only for the brute-force checks.",
         "",
         "## Run parameters (a(50))", "",
         "* Box v <= 50, f <= 96 (`combine_full.py 50 96 10`).",
         f"* Sensed counts S: `sensed_p.py 146 i 50 96`, i = 0..9, primes {', '.join(map(str, PRIMES[:10]))} "
         "(CRT in `sensed_146.json`).",
         f"* W_ref, W_rr: `run_wref.py 50 96 i`, `run_wrr.py 50 96 i`, i = 0..9, primes "
         f"{', '.join(map(str, PRIMES24[:10]))}.",
         "* CRT range: `crt_bound.py` proves 2 R(v,f) < P for both tables, R = rooted c-nets (README §7).",
         "",
         "## OEIS b-files in `data/`", "",
         "Copies of the OEIS b-files the checks read, as downloaded in October 2026 (OEIS content, CC BY-SA 4.0).",
         "`b000944.txt` is the OEIS file before the extension (a(1)-a(18)), not the new `b000944.txt` at the",
         "top of this folder.", "",
         "| file | source | terms | SHA-256 |", "|---|---|---|---|"]
    for name in DATA:
        p = HERE / "data" / name
        a = "A" + name[1:7]
        L.append(f"| `data/{name}` | https://oeis.org/{a}/{name} | {lines(p)} | `{sha(p)}` |")
    L += ["", "## Residues of the a(50) run (`formula/`)", "",
          "`combine_full.py 50 96 10` reads these and `sensed_146.json`; it rebuilds `Ufull_50_96.json` "
          "(published as `results/U_50_96.json`).", "",
          "| file | SHA-256 |", "|---|---|"]
    for pat in ("wref_50_96_*.npy", "wrr_50_96_*.npy", "sensed_146_*.npy"):
        for p in sorted(F.glob(pat)):
            L.append(f"| `formula/{p.name}` | `{sha(p)}` |")
    L += ["", "## Outputs", "", "| file | SHA-256 |", "|---|---|"]
    for rel in ("b000944.txt", "results/U_50_96.json", "results/sensed_146.json", "results/a109_80.json"):
        L.append(f"| `{rel}` | `{sha(HERE / rel)}` |")
    (HERE / "MANIFEST.md").write_text("\n".join(L) + "\n", encoding="utf-8", newline="\n")
    print("wrote MANIFEST.md")


if __name__ == "__main__":
    main()
