#!/bin/bash
# Main runs plus reflected cross-checks (reflection in the x axis maps
# Q3->Q1 to Q2->Q4, Q3->Q4 to Q2->Q1, Q3->Q3 to Q2->Q2).
N=${1:-36}
B=rust/target/release/meander2.exe
mkdir -p results
$B $N q1            > results/q1_$N.txt
$B $((N+1)) q4      > results/q4_$((N+1)).txt
$B $N q3            > results/q3_$N.txt
$B $N q4 start=q2   > results/q2q4_$N.txt
$B $((N+1)) q1 start=q2 > results/q2q1_$((N+1)).txt
$B $N q2 start=q2   > results/q2q2_$N.txt
