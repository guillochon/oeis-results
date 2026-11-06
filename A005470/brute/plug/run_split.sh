#!/bin/bash
# Usage: run_split.sh n parts outfile   (run from WSL; plantri_zk in ~/src/plantri_dl/plantri58)
n=$1; k=$2; out=$3
cd ~/src/plantri_dl/plantri58
for r in $(seq 0 $((k-1))); do nice -n 10 ./plantri_zk -pG $n $r/$k 2>/dev/null > /tmp/zk_${n}_$r.txt & done
wait
cat /tmp/zk_${n}_*.txt > "$out"
