#!/bin/bash

echo "Cleaning folder form output of previous runs of IMEX_SfloW2D"

file="TEST_CF*"

rm -f $file

file="topo*"

rm -f $file

file="*bak"

rm -f $file

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

