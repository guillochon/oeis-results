#!/bin/bash
# Re-run a small rep range several times under a timeout (hang diagnosis).
#   retest.sh FROM TO [REPEATS]
cd "$(dirname "$0")/.."
for i in $(seq 1 ${3:-3}); do
  rm -f /tmp/h.txt
  timeout 60 ./cuda/a8 rust/results/ylo_reps_d6.bin --from $1 --to $2 --batch 4 --sym --out /tmp/h.txt > /tmp/h.log 2>&1
  echo "try $i: exit=$? $(tail -1 /tmp/h.log)"
done
