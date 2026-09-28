#!/bin/bash
# GPU utilisation as reported inside WSL while a short normal run is going.
cd "$(dirname "$0")/.."
rm -f /tmp/u.txt
./cuda/a8 rust/results/ylo_reps_d6.bin --from 700000 --to 700100 --batch 4 --sym --out /tmp/u.txt > /dev/null 2>&1 &
sleep 3
for i in 1 2 3 4; do nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader; sleep 1; done
wait
