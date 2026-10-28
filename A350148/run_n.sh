#!/usr/bin/env bash
# Count every zipping mode for one grid size, one mode at a time at below-normal priority (resumable:
# modes already in results/nN.txt are skipped).   usage: bash run_n.sh N [big [THREADS [MEM_GB]]]
# With "big", uses the sharded multithreaded engine (rust/src/bin/big.rs) built in target-dev.
set -u
cd "$(dirname "$0")"
N=$1
OUT=results/n$N.txt
python make_bounds.py "$N" > /dev/null
touch "$OUT"
for f in bounds/n${N}_*.txt; do
    grep -q "^$f " "$OUT" && continue
    if [ "${2:-}" = big ]; then
        # row-by-row progress goes to results/nN.rows; the result line has the same format as `rust count`
        res=$(cmd //c "start /b /belownormal /wait rust\\target-dev\\release\\big.exe ${f//\//\\} --threads ${3:-6} --mem-gb ${4:-6}" 2>> "results/n$N.rows" | tr '\n' ' ')
    else
        res=$(cmd //c "start /b /belownormal /wait rust\\target\\release\\rust.exe count ${f//\//\\}" 2>&1 | tr '\n' ' ')
    fi
    [ -n "$res" ] || { echo "no result for $f" >&2; exit 1; }
    echo "$f $res" >> "$OUT"
    echo "$f $res"
done
echo "done n=$N"
