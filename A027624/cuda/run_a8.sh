#!/bin/bash
# Full a(8) run on the GPU (inside WSL). Resumable: rerun the same command after an interruption.
#   bash run_a8.sh [--nosym] [--fg]
# From PowerShell, use --fg (the job then lives as long as that wsl.exe session, and WSL does not
# idle the distro away under it):
#   wsl -d Ubuntu-22.04 -- bash /mnt/c/Users/guill/oeis-sniper/attacks/A027624/cuda/run_a8.sh --fg
# Without --fg the job is started in the background with nohup (for use from a WSL shell).
# Progress: attacks/A027624/logs/a8_gpu.log     Result: python attacks/A027624/cuda/aggregate.py
#
# The orbits are processed in chunks. The GPU program exits by itself if a batch stalls (the GPU
# occasionally stops completing work under WSL; exit code 3), and each chunk also has a timeout;
# either way the chunk is rerun and resumes (finished orbits are kept in the output file).
cd "$(dirname "$0")/.."
mkdir -p results logs
SYM=--sym; OUT=results/a8_sym.txt; REF=logs/cpu6_reference_sym.txt; FG=0
for arg in "$@"; do
  case "$arg" in
    --nosym) SYM=; OUT=results/a8_plain.txt; REF=logs/cpu6_reference.txt ;;
    --fg) FG=1 ;;
    *) echo "unknown option $arg" >&2; exit 1 ;;
  esac
done
if pgrep -f "cuda/a8 rust/results" > /dev/null; then echo "an a(8) job is already running (see status.sh)" >&2; exit 1; fi
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -o cuda/a8 cuda/a8.cu || exit 1
# guard: the build must reproduce the CPU reference sums before a long run
./cuda/a8 rust/results/ylo_reps_d6.bin --check $REF --batch 4 $SYM > /dev/null 2>&1 \
  || { echo "reference check failed; not starting" >&2; exit 1; }

N=1228158; CH=4000
run_chunks() {
  for ((s = 0; s < N; s += CH)); do
    for try in $(seq 1 50); do
      timeout 3600 ./cuda/a8 rust/results/ylo_reps_d6.bin --from $s --to $((s + CH)) --batch 4 $SYM --out $OUT && break
      echo "chunk $s: attempt $try failed or timed out; retrying"
    done
    echo "[$(date '+%F %T')] reps done through $((s + CH)) of $N"
  done
  echo "[$(date '+%F %T')] finished"
}

if [ $FG = 1 ]; then
  echo "running in the foreground; output $OUT; log attacks/A027624/logs/a8_gpu.log (Ctrl+C to stop, rerun to resume)"
  run_chunks >> logs/a8_gpu.log 2>&1
else
  export -f run_chunks; export N CH SYM OUT
  nohup bash -c run_chunks >> logs/a8_gpu.log 2>&1 &
  echo "started (pid $!); output $OUT; log: attacks/A027624/logs/a8_gpu.log"
fi
