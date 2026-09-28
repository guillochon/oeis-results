#!/bin/bash
# gprun.sh "gp commands"   (runs in the attack folder)
cd "$(dirname "$0")/.."
echo "$1" | gp -q --default parisizemax=4000000000
