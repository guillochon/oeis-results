#!/bin/bash
# Build, then check the GPU against both CPU references (plain and swap-weighted), then time it.
cd "$(dirname "$0")"
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -Xptxas -v -o a8 a8.cu 2>&1 | grep -o "Used [0-9]* registers" | head -1
./a8 ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference.txt --batch 4 2>&1 | tail -1
[ -f ../logs/cpu6_reference_sym.txt ] && ./a8 ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference_sym.txt --batch 4 --sym 2>&1 | grep -v "reps loaded"
for f in 0 700000 1228000; do
  rm -f /tmp/t.txt
  echo -n "sym bench reps $f..: "; ./a8 ../rust/results/ylo_reps_d6.bin --from $f --to $((f+64)) --batch 4 --sym --out /tmp/t.txt 2>&1 | tail -1
done
