#!/bin/bash
# Repeat the reference check N times with a timeout, to catch intermittent hangs.
cd "$(dirname "$0")"
N=${1:-15}; ok=0; bad=0; hang=0
for i in $(seq 1 $N); do
  r=$(timeout 60 ./a8 ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference_sym.txt --batch 4 --sym 2>&1 | tail -1)
  rc=$?
  if echo "$r" | grep -q "ALL MATCH"; then ok=$((ok+1)); elif [ -z "$r" ]; then hang=$((hang+1)); echo "run $i: hang/timeout"; else bad=$((bad+1)); echo "run $i: $r"; fi
done
echo "stress: $ok ok, $bad mismatched, $hang hung (of $N)"
