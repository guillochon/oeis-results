#!/bin/bash
cd "$(dirname "$0")/.."
M=${1:-2000000}
time bash tools/run_h3.sh $M 1000000 2
cat data/h3/*.done
echo "M=$M; read(\"tools/h3_validate.gp\");" | gp -q --default parisizemax=4000000000
