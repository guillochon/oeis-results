#!/bin/bash
# stage6_par.sh NAME SIGN
cd "$(dirname "$0")/.."
NAME=$1; SIGN=$2
for f in data/s6/${NAME}_??.txt; do
  printf 'infile="%s"; sgn=%s; outfile="%s.out";\nread("tools/stage6_run.gp");\n' $f $SIGN $f > $f.gp
  gp -q --default parisizemax=4000000000 $f.gp < /dev/null > $f.log 2>&1 &
done
wait
cat data/s6/${NAME}_??.txt.out > data/stage6_${NAME}.txt
cat data/s6/${NAME}_??.txt.log | grep -v Warning
