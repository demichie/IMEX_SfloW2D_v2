#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of SW_VAR_DENS_MODEL"

for file in exampleBOX*
do
    if [ -f $file ] ; then
        rm $file
    fi
done

file="SW_VAR_DENS_MODEL.inp"

if [ -f $file ] ; then
    rm $file
fi

file="example_BOX_0000.q_2d"

if [ -f $file ] ; then
    rm $file
fi

file="topography_dem.asc"

if [ -f $file ] ; then
    rm $file
fi


