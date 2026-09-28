#!/bin/bash
# runfile.sh file.gp logfile : run a gp script in the attack folder, log stdout+stderr
cd "$(dirname "$0")/.."
gp -q --default parisizemax=4000000000 "$1" < /dev/null > "$2" 2>&1
