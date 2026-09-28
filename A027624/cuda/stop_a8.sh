#!/bin/bash
# Stop a running a(8) job (the chunk loop and the GPU program).
pkill -f "run_a8.sh" ; pkill -f run_chunks ; pkill -f "cuda/a8 rust/results" ; sleep 1
ps -eo pid,cmd | grep -E "cuda/a8 rust|run_a8.sh|run_chunks" | grep -v grep || echo "a(8) job stopped"
