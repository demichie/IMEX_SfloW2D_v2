#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of SW_VAR_DENS_MODEL"

file="TEST_BW*"

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

