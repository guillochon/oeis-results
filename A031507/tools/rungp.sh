#!/bin/bash
# usage: rungp.sh file.gp  (runs from the attack folder)
cd "$(dirname "$0")/.."
gp -q --default parisizemax=4000000000 "$@" < /dev/null
