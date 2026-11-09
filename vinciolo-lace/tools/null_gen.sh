#!/bin/bash
# 20 density-matched null designs per panel: interior size n-2, padded with a frame later.
mkdir -p null
python3 - <<'PY' > null/jobs.txt
import json
for p in json.load(open('panels.json')):
    I=[r[1:-1] for r in p['rows'][1:-1]]; d=sum(c=='1' for r in I for c in r)/len(I)**2
    print(p['name'], p['n']-2, round(d,4))
PY
while read name m d; do ./tools/null_sampler $m $d 7 20 > null/$name.txt 2>/dev/null; done < null/jobs.txt
