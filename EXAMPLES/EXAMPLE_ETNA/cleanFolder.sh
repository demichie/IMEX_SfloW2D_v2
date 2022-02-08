#!/bin/bash
clear
echo "Cleaning folder form output of previous runs of SW_VAR_DENS_MODEL"

file="exampleEtna*"

rm -f $file

file="dem_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_interfaces_esri.asc"

if [ -f $file ] ; then
    rm $file
fi

file="dem_esri_nodata.asc"

if [ -f $file ] ; then
    rm $file
fi

file="./DEM/ellipsoid.asc"

if [ -f $file ] ; then
    rm $file
fi

file="./DEM/init_volume.asc"

if [ -f $file ] ; then
    rm $file
fi

file="./DEM/new_dem.asc"

if [ -f $file ] ; then
    rm $file
fi

folder="__pycache__"

if [ -d $folder ] ; then
    rm -rf $file
fi
