#!/bin/bash
# Build fixcount and write logs/fixcount_<n>[_<mode>].txt (WSL or any Linux with gcc).
#   bash run_fixcount.sh all 2 3 4 5 6 7 8 9     (every pair of cycle types)
#   bash run_fixcount.sh fpf 10 11               (fixed-point-free pairs only)
#   SKIP2=skip2 bash run_fixcount.sh fpf 10      (... and not the (2^5, 2^5) pair)
cd "$(dirname "$0")"
gcc -O2 -o fixcount fixcount.c || exit 1
mode=$1; shift
for n in "$@"; do
  suffix=""; [ "$mode" = all ] || suffix="_$mode"
  /usr/bin/time -f "n=$n $mode: %e s" ./fixcount "$n" "$mode" $SKIP2 > "logs/fixcount_${n}${suffix}.txt"
done
