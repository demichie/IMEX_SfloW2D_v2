# Depth-averaged gas-particles model

[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_gold.png)](https://api.eu.badgr.io/public/assertions/VYyqk62SQfq4Hx1E3CsHHg "SQAaaS gold badge achieved")

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-gold-yellow)](https://api.eu.badgr.io/public/assertions/VYyqk62SQfq4Hx1E3CsHHg "SQAaaS gold badge achieved")

[![DOI](https://zenodo.org/badge/218571198.svg)](https://zenodo.org/doi/10.5281/zenodo.7476736)

Shallow water model for multiphase flow (gas+particles) with density of gas
temperature-dependent.

To compile the code you need a Fortran compiler and the NetCDF library for Fortran. You can install both with anaconda, by creating an anaconda environment and activating it:

> conda create -n fortran_env conda-forge::gfortran_linux-64  sysroot_linux-64  make  conda-forge::netcdf-fortran  conda-forge::liblapack  conda-forge::libblas
> 
> conda activate fortran_env

To compile:

> autoreconf -i

Then on linux (replace USERNAME with you user account):

> ./configure --with-netcdf=/home/USERNAME/anaconda3/envs/fortran_env

On OSX (replace USERNAME with you user account):

> ./configure --with-netcdf=/USERS/USERNAME/anaconda3/envs/fortran_env

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: -fopenmp
2) with intel: -qopenmp

> make

> make install

The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.

## Docker container

If you do not have a compiler on your system, there is a Docker container
with the executable of the latest version of the model at the following
link:

<https://github.com/demichie/IMEX_SfloW2D_v2/pkgs/container/imex_sflow2d_v2>

If you have docker installed on your computer, you can download the container
from the commant line with:

> docker pull ghcr.io/demichie/imex_sflow2d_v2:sha256-080bf69dd96c57482ca27fa37096e96b9b15813ca5d7a8a5736ea9a0599884f9

### LINUX/MAC

Create a folder for your simulation with all the input files and then run the
container with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those
files to have .png files with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those
files to have a netCDF4 file with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf

### WINDOWS COMMAND LINE

Create a folder for your simulation with all the input files and then run
the container with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those
files to have .png files with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those
files to have a netCDF4 file with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf

### WINDOWS POWERSHELL

Create a folder for your simulation with all the input files and then run
the container with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those
files to have .png files with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those
files to have a netCDF4 file with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf

