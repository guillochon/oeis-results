#!/bin/bash
# Usage: run_gpu.sh N NPRIMES   -- formula rows on the GPU (KF_GPU=2, two streams), then the tower.
cd "$(dirname "$0")"
N=$1; NP=$2
export KF_GPU=2 PYTHONWARNINGS=ignore
stream() { for i in "$@"; do [ -s "kf_${N}_done_$i" ] && continue; cmd //c "start /b /wait /belownormal python kformula.py $N $i" > kf${N}_$i.log 2>&1 && echo ok > "kf_${N}_done_$i" || return 1; done; }
evens=$(seq 0 2 $((NP-1))); odds=$(seq 1 2 $((NP-1)))
stream $evens & A=$!
stream $odds & B=$!
wait $A $B
grep -l Traceback kf${N}_*.log && { grep -A8 Traceback kf${N}_*.log; exit 1; }
rm -f kf_${N}_done_*
cmd //c "start /b /wait /belownormal python tower_gpu.py $N $NP" > tower_$N.log 2>&1
tail -30 tower_$N.log
