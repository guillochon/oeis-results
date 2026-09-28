#!/bin/bash
# Status of the a(8) run: processes, progress, GPU.
cd "$(dirname "$0")/.."
ps -eo pid,etime,pcpu,cmd | grep -E "cuda/a8 rust|run_a8.sh|run_chunks|timeout [0-9]" | grep -v grep | cut -c1-150
echo "orbits done: $(cat results/a8_sym.txt 2>/dev/null | wc -l) / 1228158"
tail -3 logs/a8_gpu.log
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader
