#!/bin/bash
# Run unil for primes FIRST..LAST of primes.txt, JOBS at a time, at low priority:  bash run.sh N JOBS [FIRST [LAST]]
# (WSL: gcc -O2 -march=native -o unil unil.c)
cd "$(dirname "$0")"
N=$1; JOBS=${2:-4}; FIRST=${3:-0}; LAST=${4:-999999}
mkdir -p results
nl -v0 primes.txt | sed -n "$((FIRST + 1)),$((LAST + 1))p" | xargs -P "$JOBS" -L1 bash -c 'nice -n 19 ./unil '"$N"' $1 > results/n'"$N"'_p$0.txt 2> results/n'"$N"'_p$0.log'
