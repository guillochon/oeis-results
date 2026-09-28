#!/bin/bash
# Tests for the host-side robustness: build, reference check, resume from a truncated line, and
# the stall watchdog (with --stall 0.001 every batch "stalls").
cd "$(dirname "$0")/.."
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -o cuda/a8 cuda/a8.cu || exit 1
./cuda/a8 rust/results/ylo_reps_d6.bin --check logs/cpu6_reference_sym.txt --batch 4 --sym 2>&1 | tail -1
# resume: 8 reps, then truncate the last line and rerun: only that orbit should be redone
rm -f /tmp/r.txt
./cuda/a8 rust/results/ylo_reps_d6.bin --from 1668 --to 1676 --batch 4 --sym --out /tmp/r.txt > /dev/null 2>&1
cp /tmp/r.txt /tmp/r_full.txt
head -c -10 /tmp/r_full.txt > /tmp/r.txt      # cut the last line short (no newline)
./cuda/a8 rust/results/ylo_reps_d6.bin --from 1668 --to 1676 --batch 4 --sym --out /tmp/r.txt 2>&1 | grep -E "reps to do|dropping"
cmp <(sort /tmp/r.txt) <(sort /tmp/r_full.txt) && echo "resume after truncation: file identical to the uninterrupted run"
# watchdog
timeout 60 ./cuda/a8 rust/results/ylo_reps_d6.bin --from 700000 --to 700004 --batch 4 --sym --stall 0.001 --out /tmp/w.txt 2>&1 | grep STALL
echo "exit code with forced stall: ${PIPESTATUS[0]}"
