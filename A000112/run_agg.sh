#!/bin/bash
# Run in WSL: bash run_agg.sh M PARTS
# Aggregates of the M-point posets for the block engine, in PARTS parallel genposetg slices
# (genposetg's "m x y" option); then concatenates the parts into aggs/aggM.txt.
set -e
cd "$(dirname "$0")"
G=~/src/nauty2_8_9/genposetg
M=$1; PARTS=${2:-6}
mkdir -p aggs logs
for x in $(seq 0 $((PARTS-1))); do
  ( "$G" "$M" t m "$x" "$PARTS" 2>/dev/null | rust/target/release/a000112_run.exe agg "$M" > "aggs/agg${M}_p${x}.txt" 2> "logs/agg${M}_p${x}.log" ) &
done
wait
cat aggs/agg${M}_p*.txt > aggs/agg${M}.txt
grep -h "posets" logs/agg${M}_p*.log
