#!/bin/bash
cd "$(dirname "$0")/.."
echo "K=$1; read(\"tools/cp_dump.gp\");" | gp -q --default parisizemax=4000000000
