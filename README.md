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

or, on a osx computer:

> conda create -n fortran_env conda-forge::gfortran make conda-forge::netcdf-fortran conda-forge::liblapack conda-forge::libblas
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

