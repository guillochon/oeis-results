#!/bin/bash
# Compare kernel variants: correctness (swap-weighted reference) and speed on a uniform sample.
cd "$(dirname "$0")"
for v in "" "-DTHREADS=512" "-DTHREADS=384"; do
  /usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w $v -Xptxas -v -o /tmp/a8v a8.cu 2>&1 | grep -o "Used [0-9]* registers.*smem" | head -1
  c=$(/tmp/a8v ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference_sym.txt --batch 4 --sym 2>&1 | tail -1)
  rm -f /tmp/t.txt
  r=$(/tmp/a8v ../rust/results/ylo_reps_d6.bin --stride 4096 --batch 4 --sym --out /tmp/t.txt 2>&1 | tail -1)
  printf "%-16s check [%s]  sample [%s]\n" "[$v]" "$c" "$r"
done
