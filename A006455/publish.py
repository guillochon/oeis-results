"""Publish A006455's shareable results to github.com/guillochon/oeis-results/tree/main/A006455.

    python attacks/A006455/publish.py --dry-run   # show what would change
    python attacks/A006455/publish.py             # commit and push
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from publish import cli, lean_project  # noqa: E402

FILES = [
    "README.md",
    "layman.md",
    "b006455.txt",
    "*.py",
    "gen_posets.sh",
    "rust/Cargo.toml",
    "rust/Cargo.lock",
    "rust/src/*.rs",
    "logs/*.log",
]

# Sequences these files may mention (poset counts and related entries in the OEIS).
ALLOWED_REFS = {"A000112", "A001035", "A121337", "A323502", "A000798"}

if __name__ == "__main__":
    cli("A006455", HERE, FILES, ALLOWED_REFS, extra=lean_project(["A006455"]))
