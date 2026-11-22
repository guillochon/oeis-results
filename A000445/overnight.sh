#!/bin/bash
# Overnight a(9) ladder: portfolio search (THREADS threads, default 4) at increasing N until STOP
# (default 06:55 tomorrow). Each SAT witness is checked by verify_witness_chunked.py and becomes the hint
# for the next rung; a TIMEOUT or UNSAT ends the ladder.
# usage: [THREADS=t] [PREV=hint_witness] [STOP=epoch] overnight.sh [N1 N2 ...]  (THREADS=1 runs only the hinted search)
cd "$(dirname "$0")"
B=rust/target/release/cpr.exe
T=${THREADS:-4}
STOP=${STOP:-$(date -d "tomorrow 06:55" +%s)}
prev=${PREV:-results/witness_k9_3e8.txt}
LADDER=${*:-350000000 400000000 450000000 500000000}
out=results/overnight_k9.log
echo "$(date +%H:%M) start: threads=$T ladder=$LADDER stop=$(date -d @$STOP +%H:%M)" >> $out
for N in $LADDER; do
  left=$((STOP - $(date +%s)))
  [ $left -lt 600 ] && { echo "$(date +%H:%M) out of time before N=$N" >> $out; break; }
  log=results/portfolio_k9_$N.log
  cmd //c start //b //wait //belownormal "$B" portfolio 9 $N $left results/witness_k9_$N.txt $T $prev > $log 2>&1
  res=$(grep -aE "^(SAT|UNSAT|TIMEOUT)" $log)
  echo "$(date +%H:%M) N=$N -> $res | $(grep -aE '^k=' $log | tail -1)" >> $out
  case "$res" in
    SAT*)
      cmd //c start //b //wait //belownormal python verify_witness_chunked.py 9 $N results/witness_k9_$N.txt \
        --extend $((N + 20000000)) > results/verify_k9_$N.log 2>&1
      echo "  verify: $(tail -2 results/verify_k9_$N.log | tr '\n' ' ')" >> $out
      prev=results/witness_k9_$N.txt ;;
    *) break ;;
  esac
done
echo "$(date +%H:%M) done" >> $out
