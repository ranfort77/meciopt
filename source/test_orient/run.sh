#!/bin/bash

wfu="/home/user/ahs/meciopt/orient_tests"
scdir="/work/tmp/ahs/molpro_$$"
prefix=`printf $1 | cut -d '.' -f 1`
out=${prefix}.out

molpro -s --no-xml-output -d $scdir -W $wfu -o $out $1
