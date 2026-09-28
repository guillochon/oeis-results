#!/bin/bash
# GPU throughput on a few blocks of reps spread over the rep list.
cd "$(dirname "$0")/.."
for f in 0 300000 700000 1100000 1228000; do
  rm -f /tmp/bench_$f.txt
  echo -n "reps $f..$((f+64)): "
  ./cuda/a8 rust/results/ylo_reps_d6.bin --from $f --to $((f+64)) --batch ${1:-4} --out /tmp/bench_$f.txt 2>&1 | tail -1
done
