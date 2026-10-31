#!/bin/bash
# usage: bash run.sh N  -> results/nN_p{0..4}.txt, one prime after another, 6 threads
N=$1
for p in 0 1 2 3 4; do
  ./rust/mtg_run.exe $N $p 6 > results/n${N}_p$p.tmp 2> results/n${N}_p$p.log && mv results/n${N}_p$p.tmp results/n${N}_p$p.txt
done
