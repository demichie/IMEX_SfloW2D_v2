#!/bin/bash
cd /home/user_sw/SW_RUNS

case $1 in
   plot_overlay)
      echo "plot_overlay"
      touch matplotlibrc
      echo "backend : agg" >> matplotlibrc      
      python3 /home/user_sw/IMEX_SfloW2D-master/UTILS/plot_overlay.py
      rm matplotlibrc;;
   p2d_to_netcdf)
      echo "netcdf"
      python3 /home/user_sw/IMEX_SfloW2D-master/UTILS/p2d_to_netCDF4.py;;
   run)   
      echo "run"
      /home/user_sw/IMEX_SfloW2D-master/bin/IMEX_SfloW2D;;
   *)
      echo "Please select: run, plot_overlay, p2d_to_netcdf"
      echo"";;
esac

#
# pwd

