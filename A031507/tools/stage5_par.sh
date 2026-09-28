#!/bin/bash
# stage4_par.sh NAME SIGN PMAX
cd "$(dirname "$0")/.."
NAME=$1; SIGN=$2; PMAX=${3:-1000}
for f in data/s5/${NAME}_??.txt; do
  printf 'infile="%s"; sgn=%s; outfile="%s.out"; PMAX=%s;\nread("tools/stage5_run.gp");\n' $f $SIGN $f $PMAX ${TGT:-6} > $f.gp
  gp -q --default parisizemax=4000000000 $f.gp < /dev/null > $f.log 2>&1 &
done
wait
cat data/s5/${NAME}_??.txt.out > data/stage5_${NAME}.txt
cat data/s5/${NAME}_??.txt.log | grep -v Warning
