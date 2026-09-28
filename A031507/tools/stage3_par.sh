#!/bin/bash
# stage3_par.sh NAME SIGN : run stage 3 on data/s2/NAME_part_??.out in parallel
cd "$(dirname "$0")/.."
NAME=$1; SIGN=$2
for f in data/s2/${NAME}_part_??.out; do
  printf 'infile="%s"; sgn=%s; outfile="%s.s3";\nread("tools/stage3_run.gp");\n' $f $SIGN $f > $f.s3.gp
  gp -q --default parisizemax=4000000000 $f.s3.gp < /dev/null > $f.s3.log 2>&1 &
done
wait
cat data/s2/${NAME}_part_??.out.s3 > data/stage3_${NAME}.txt
cat data/s2/${NAME}_part_??.out.s3.log
