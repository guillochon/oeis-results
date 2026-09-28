#!/bin/bash
# one_chunk.sh lo hi sizemax
cd "$(dirname "$0")/.."
out=data/h3/h3_$(printf %010d $1).txt
rm -f $out $out.done
/usr/bin/time -v bash -c "echo 'lo=$1; hi=$2; out=\"$out\"; read(\"tools/h3chunk.gp\");' | gp -q --default parisizemax=$3" 2>&1 | grep -E "Maximum resident|Elapsed|overflow"
cat $out.done
