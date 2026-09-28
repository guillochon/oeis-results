#!/bin/bash
# Nsight Compute profile of one kernel launch (1 orbit, all 65536 slices), with source-line info.
# (WSL: no PM sampling, so explicit sections instead of --set full.)
cd "$(dirname "$0")"
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -lineinfo -o /tmp/a8p a8.cu || exit 1
rm -f /tmp/p.txt
SECTIONS="${SECTIONS:---section SpeedOfLight --section Occupancy --section LaunchStats --section WarpStateStats --section InstructionStats --section ComputeWorkloadAnalysis --section MemoryWorkloadAnalysis --section SchedulerStats --section SourceCounters}"
/usr/local/cuda-13.3/bin/ncu $SECTIONS --import-source yes -k slice_kernel -c 1 -f -o /tmp/a8prof \
  /tmp/a8p ../rust/results/ylo_reps_d6.bin --from 700000 --to 700001 --batch 1 --sym --out /tmp/p.txt > /tmp/ncu.log 2>&1
grep -E "ERROR|PROF" /tmp/ncu.log | head -5
/usr/local/cuda-13.3/bin/ncu --import /tmp/a8prof.ncu-rep --page details 2>&1 | grep -v "^ *-*$" | head -${LINES_OUT:-220}
