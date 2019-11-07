#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of SW_VAR_DENS_MODEL"

file="example2D*"

rm -f $file

file="example_2D_0000.q_2d"

if [ -f $file ] ; then
    rm $file
fi


file="SW_VAR_DENS_MODEL.inp"

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

