#!/bin/bash

sed 's/^>/@/'  ${@:-/dev/stdin} | paste - - | mawk -vFS="\t" -vOFS="\n" '{print $1,$2,"+";gsub(".","a",$2);print $2}'

