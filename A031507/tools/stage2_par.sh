#!/bin/bash
# stage2_par.sh NAME SIGN P : split data/cand_NAME_7_56877643.txt into P parts and run stage 2 in parallel
cd "$(dirname "$0")/.."
NAME=$1; SIGN=$2; P=${3:-10}
IN=data/cand_${NAME}_7_56877643.txt
mkdir -p data/s2
split -n l/$P -d -a 2 $IN data/s2/${NAME}_part_
for f in data/s2/${NAME}_part_??; do
  printf 'infile="%s"; sgn=%s; outfile="%s.out";\nread("tools/stage2_run.gp");\n' $f $SIGN $f > $f.gp
  gp -q --default parisizemax=4000000000 $f.gp < /dev/null > $f.log 2>&1 &
done
wait
cat data/s2/${NAME}_part_??.out > data/stage2_${NAME}.txt
cat data/s2/${NAME}_part_??.log
