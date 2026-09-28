#!/bin/bash
for c in mwrank magma sage; do command -v $c || echo "no $c"; done
dpkg -l 2>/dev/null | grep -iE "eclib|mwrank|pari" | awk '{print $2, $3}'
