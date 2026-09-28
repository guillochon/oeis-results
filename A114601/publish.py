"""Publish A114601's shareable results to github.com/guillochon/oeis-results/tree/main/A114601.

    python attacks/A114601/publish.py --dry-run   # show what would change
    python attacks/A114601/publish.py             # commit and push
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from publish import cli, lean_project  # noqa: E402

FILES = [
    "README.md",
    "layman.md",
    "proofs.md",
    "b114601.txt",
    "*.py",
    "rust/Cargo.toml",
    "rust/Cargo.lock",
    "rust/.cargo/config.toml",
    "rust/src/*.rs",
    "logs/*.log",
]

# Sequences these files may mention (related matrix / tree / reflectable-basis counts in the OEIS).
ALLOWED_REFS = {"A086215", "A085657", "A057500", "A000272", "A000169", "A000435", "A320064", "A000616"}

if __name__ == "__main__":
    cli("A114601", HERE, FILES, ALLOWED_REFS, extra=lean_project(["A114601"]))
