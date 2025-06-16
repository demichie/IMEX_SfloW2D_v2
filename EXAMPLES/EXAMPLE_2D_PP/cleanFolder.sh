#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of IMEX_SfloW2D"

file="example2D*"

rm -f $file

file="example_2D_0000.q_2d"

if [ -f $file ] ; then
    rm $file
fi


file="IMEX_SfloW2D.inp"

if [ -f $file ] ; then
    rm $file
fi

file="topography_dem.asc"

if [ -f $file ] ; then
    rm $file
fi

file="pile.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

