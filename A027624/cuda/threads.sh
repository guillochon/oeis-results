#!/bin/bash
# Block-size experiment: correctness check + timing for several THREADS values.
cd "$(dirname "$0")"
for t in ${THREADS_LIST:-256 512}; do
  /usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -DTHREADS=$t -Xptxas -v -o /tmp/a8t a8.cu 2>&1 | grep -o "Used [0-9]* registers" | head -1
  c=$(/tmp/a8t ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference.txt --batch 4 2>&1 | tail -1)
  rm -f /tmp/t.txt
  r=$(/tmp/a8t ../rust/results/ylo_reps_d6.bin --from 700000 --to 700064 --batch 4 --out /tmp/t.txt 2>&1 | tail -1)
  echo "THREADS=$t: check [$c]  bench [$r]"
done
