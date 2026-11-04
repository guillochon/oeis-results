#!/bin/bash
# Run in WSL: bash run_stat2.sh R "L1,L2/L1,L2,L3/..." PARTS TAG ["x1 x2 ..."]
# Raw Burnside sums for several interacting-length lists over the R-point posets, in PARTS
# parallel genposetg slices (all of them, or only the listed slice indices); combine with
# combine_stats.py.
set -e
cd "$(dirname "$0")"
G=~/src/nauty2_8_9/genposetg
R=$1; KEYS=$2; PARTS=${3:-6}; TAG=${4:-run}; SLICES=${5:-$(seq 0 $((PARTS-1)))}
mkdir -p stats_raw logs
for x in $SLICES; do
  ( "$G" "$R" t m "$x" "$PARTS" 2>/dev/null | A000112_CKPT="stats_raw/stat${R}_${TAG}_p${x}.ckpt" WSLENV=A000112_CKPT/p rust/target/release/${EXE:-a000112_s3}.exe stat "$R" "$KEYS" \
      > "stats_raw/stat${R}_${TAG}_p${x}.txt" 2> "logs/stat${R}_${TAG}_p${x}.log" ) &
done
wait
grep -h "posets" logs/stat${R}_${TAG}_p*.log
