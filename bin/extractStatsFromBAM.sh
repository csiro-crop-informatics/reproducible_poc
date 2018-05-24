#!/bin/bash

awk -vOFS="\t" '{split($1,sim,"|");if(sim[4]==$3 && sim[5]==$4-1){count++}};END{print count,NR,count/NR}' ${@:-/dev/stdin}


