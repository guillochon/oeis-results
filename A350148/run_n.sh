#!/usr/bin/env bash
# Count every zipping mode for one grid size, one mode at a time at below-normal priority (resumable:
# modes already in results/nN.txt are skipped).   usage: bash run_n.sh N
set -u
cd "$(dirname "$0")"
N=$1
OUT=results/n$N.txt
python make_bounds.py "$N" > /dev/null
touch "$OUT"
for f in bounds/n${N}_*.txt; do
    grep -q "^$f " "$OUT" && continue
    res=$(cmd //c "start /b /belownormal /wait rust\\target\\release\\rust.exe count ${f//\//\\}" 2>&1 | tr '\n' ' ')
    echo "$f $res" >> "$OUT"
    echo "$f $res"
done
echo "done n=$N"
