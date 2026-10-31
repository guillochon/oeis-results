#!/bin/bash
# usage: bash run2.sh N  -> results/d{N}_p{0..4}.txt ("n A000080 A048781" mod p), after any running mtg job
N=$1
while tasklist 2>/dev/null | grep -qi "mtg_run.exe"; do sleep 60; done
for p in 0 1 2 3 4; do
  ./rust/mtg2_run.exe $N $p 6 > results/d${N}_p$p.tmp 2> results/d${N}_p$p.log && mv results/d${N}_p$p.tmp results/d${N}_p$p.txt
done
