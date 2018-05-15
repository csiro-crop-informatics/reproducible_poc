#!/bin/bash

sed 's/^>/@/'  ${@:-/dev/stdin} | paste - - | awk -vFS="\t" -vOFS="\n" '{print $1,$2,"+";gsub(".","a",$2);print $2}'

