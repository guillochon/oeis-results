#!/bin/bash
# Build the a(8) GPU program (run inside WSL) and check it against the CPU reference sums.
cd "$(dirname "$0")"
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -Xptxas -v -o a8 a8.cu 2>&1 | grep -E "error|registers" 
[ "$1" = "check" ] && ./a8 ../rust/results/ylo_reps_d6.bin --check ../logs/cpu6_reference.txt --batch 4
exit 0
