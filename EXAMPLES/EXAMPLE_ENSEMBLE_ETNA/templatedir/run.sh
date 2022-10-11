#!/bin/bash
echo $PWD

mv input_ellipsoid.inp input_ellipsoid.py
python create_input_ellipsoid.py > preprocessing.log

export OMP_NUM_THREADS=1
./IMEX_SfloW2D > run.log

#python plot_overlay.py > postprocessing1.log
#zip raster.zip *.asc
#rm *.asc

#python p2d_to_netCDF4.py > postprocessing2.log
#rm *.p_2d
