#!/bin/bash
# Build the unconditional h3 table for lo <= |D| <= MAX in chunks, P parallel gp processes.
# usage: run_h3.sh MAX CHUNK P [MIN]
cd "$(dirname "$0")/.."
MAX=${1:-683000000}; CH=${2:-10000000}; P=${3:-12}; MIN=${4:-1}
mkdir -p data/h3
worker() {
  local r=$1
  local lo
  for ((lo=MIN+r*CH; lo<=MAX; lo+=P*CH)); do
    local hi=$((lo+CH-1)); [ $hi -gt $MAX ] && hi=$MAX
    local out=data/h3/h3_$(printf %010d $lo).txt
    [ -f "$out.done" ] && continue
    rm -f "$out"
    echo "lo=$lo; hi=$hi; out=\"$out\"; read(\"tools/h3chunk.gp\");" | gp -q --default parisizemax=${SIZEMAX:-10000000000} >> data/h3/worker_$r.log 2>&1
  done
}
for ((r=0; r<P; r++)); do worker $r & done
wait
echo ALL DONE
