#!/bin/bash
# 2-connected planar graphs with min degree 3 on n vertices via nauty (geng | planarg).
# 3-connectivity is filtered later in Python (pickg --G3 crashes in this nauty build).
# Usage (from PowerShell): wsl -d Ubuntu-22.04 -- bash gen_poly.sh NMAX OUTDIR
N=${1:-10}
OUT=${2:-/tmp}
NAUTY=~/src/nauty2_8_9
for n in $(seq 4 "$N"); do
  $NAUTY/geng -C -d3 -q "$n" "0:$((3*n-6))" | $NAUTY/planarg -q > "$OUT/planar2c_$n.g6"
  echo "$n $(wc -l < "$OUT/planar2c_$n.g6")"
done
