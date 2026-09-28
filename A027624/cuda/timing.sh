#!/bin/bash
# Timing ablations (results are wrong): which parts of the kernel cost what.
cd "$(dirname "$0")"
NVCC=/usr/local/cuda-13.3/bin/nvcc
run() {
  $NVCC -O3 -arch=sm_86 -w $1 -o /tmp/a8v a8.cu > /dev/null 2>&1 || { echo "build failed: $1"; exit 1; }
  rm -f /tmp/t.txt
  r=$(/tmp/a8v ../rust/results/ylo_reps_d6.bin --from 700000 --to 700064 --batch 4 --out /tmp/t.txt 2>&1 | tail -1)
  printf "%-40s %s\n" "[$1]" "$r"
}
run ""
run "-DSKIP_C"
run "-DSKIP_MUL"
run "-DTERM_MIN"
run "-DTERM_MIN -DSKIP_PERA"
