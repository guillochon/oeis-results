#!/bin/bash
# Run in WSL: bash run_stat.sh R L1,L2,...   (appends the result to stats.txt)
# Streams genposetg's R-point posets into the Windows engine.
set -e
cd "$(dirname "$0")"
G=~/src/nauty2_8_9/genposetg
"$G" "$1" t 2>/dev/null | rust/target/release/a000112_run.exe stat "$1" "$2" | tee -a stats.txt
