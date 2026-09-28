#!/bin/bash
f=/usr/share/pari/doc/usersch3.tex
n=$(grep -n 'ellrank' $f | head -3); echo "$n"
l=$(grep -n 'subsec{ellrank' $f | head -1 | cut -d: -f1)
[ -z "$l" ] && l=$(grep -n 'ellrank(' $f | head -1 | cut -d: -f1)
sed -n "$((l-2)),$((l+75))p" $f
