#!/bin/bash
# Stages 2-6 for one sign, after the Rust filter (stage 1) has written the candidate file.
# usage: pipeline.sh NAME SGN CANDFILE LIM TGT [P]
#   NAME     run name (output in data/run_NAME/)
#   SGN      1 (y^2 = x^3 + k) or -1 (y^2 = x^3 - k)
#   CANDFILE lines "k sa sb d" from `a031507 filter K TGT+1`
#   LIM      only k < LIM matter (the candidate value)
#   TGT      target: prove rank <= TGT for every k < LIM
# Each stage keeps only the curves whose proven bound is still > TGT.
set -e
cd "$(dirname "$0")/.."
NAME=$1; SGN=$2; CAND=$3; LIM=$4; TGT=$5; P=${6:-10}
D=data/run_$NAME; rm -rf $D; mkdir -p $D
GP="gp -q --default parisizemax=4000000000"

par() {  # par STAGE INFILE EXTRA  -> runs tools/STAGE_run.gp on P slices, output INFILE.STAGE
  local st=$1 in=$2 extra=$3
  split -n l/$P -d -a 2 $in $in.part_
  for f in $in.part_??; do
    [ -s $f ] || { : > $f.out; continue; }
    printf 'infile="%s"; sgn=%s; outfile="%s.out"; %s\nread("tools/%s_run.gp");\n' $f $SGN $f "$extra" $st > $f.gp
    $GP $f.gp < /dev/null > $f.log 2>&1 &
  done
  wait
  cat $in.part_??.out > $in.$st
  if grep -h "error\|\*\*\*" $in.part_??.log | grep -v Warning | grep -q .; then
    echo "ERRORS in $st:"; grep -h "\*\*\*" $in.part_??.log | grep -v Warning | head; exit 1
  fi
}

awk -v L=$LIM '$1 < L' $CAND > $D/s1.txt
echo "stage 1 (3-isogeny ambient bound > $TGT, k < $LIM): $(wc -l < $D/s1.txt) curves"

par stage2 $D/s1.txt ""
awk -v T=$TGT '$6 > T {print $1, $2, $3, $4, $5}' $D/s1.txt.stage2 > $D/s2.txt
echo "stage 2 (Cassels): $(wc -l < $D/s2.txt) left"

par stage4 $D/s2.txt "PMAX=1000; TGT=$TGT;"
awk -v T=$TGT '$6 > T {print $1}' $D/s2.txt.stage4 > $D/s4_keys.txt
awk 'NR==FNR{k[$1]=1; next} ($1 in k)' $D/s4_keys.txt $D/s2.txt > $D/s4.txt
echo "stage 4 (witnesses, p <= 1000): $(wc -l < $D/s4.txt) left"

par stage5 $D/s4.txt "PMAX=0; TGT=$TGT;"
awk -v T=$TGT '$6 > T {print $1}' $D/s4.txt.stage5 > $D/s5.txt
echo "stage 5 (full-span witnesses): $(wc -l < $D/s5.txt) left"

par stage6 $D/s5.txt ""
echo "stage 6 (certified 2-descent): $(wc -l < $D/s5.txt.stage6) curves;" \
     "certification failures: $(awk '$2 != 1' $D/s5.txt.stage6 | wc -l);" \
     "R > $TGT: $(awk -v T=$TGT '$3 > T' $D/s5.txt.stage6 | wc -l)"
