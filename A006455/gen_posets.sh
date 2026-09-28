#!/bin/bash
# Unlabeled posets on n points (Hasse diagrams, digraph6, topologically labelled), via nauty genposetg.
# Run in WSL: bash gen_posets.sh 10
set -e
G=~/src/nauty2_8_9/genposetg
cd "$(dirname "$0")"; mkdir -p data
for n in $(seq 1 "${1:-9}"); do "$G" "$n" t > "data/posets$n.d6" 2>/dev/null; done
wc -l data/posets*.d6
