#!/bin/bash

echo "Cleaning folder form output of previous runs of IMEX_SfloW2D"

file="TEST_BOX*"

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

