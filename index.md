# Welcome to IMEX_SfloW2D 2.0

IMEX_SfloW2D 2.0 is a FORTRAN90 code designed to model shallow gas-particles flows over digital elevation models (DEMs) of natural terrain. The model solves for the conservation equations (mass,momentum,energy) of the mixture, and accounts for sedimentation, erosion, friction and entrainment. The system is described by an hyperbolic system of partial differential equations with relaxation and source terms. It is possible to select a simpler rheology in order to mimic the system of equations described in Kurganov and Petrova, 2007.

Several examples (1D and 2D) are provided with the package.
The code can deal with different scenarios, but its first aim is to treat gravitational flows over topographies described as digital elevation models (DEMs) in the ESRI ascii format. Moreover it is possible to save the solution as as ESRI ascii files, suitable for GIS softwares.

### Authors and Contributors

Mattia de' Michieli Vitturi (@demichie)

Tomaso Esposti Ongaro

Samantha Engwell

Brandon Keim

### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/IMEX_SfloW2D_v2/wiki](https://github.com/demichie/IMEX_SfloW2D_v2/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/IMEX_SfloW2D_v2/html/](http://demichie.github.io/IMEX_SfloW2D_v2/html/) 


### Installation and execution

Before compiling the code, please be sure that the following libraries are installed on your system:

- liblapack-dev 
- libopenblas-dev


Download the IMEX_SfloW2D_v2 package and create the executable with the following commands from a terminal:

> touch README
> 
> autoreconf
> 
> ./configure

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: -fopenmp
2) with intel: -qopenmp

> make

> make install

This will create the executable and copy it in the bin folder. You can test the executable copying it in the EXAMPLES folder and running it.

The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.

It is possible also to use Docker container with the latest version of the model:

> docker pull demichie/imex_sflow2d_v2

LINUX/MAC

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf

WINDOWS COMMAND LINE

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf

WINDOWS POWERSHELL

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/imex_sflow2d_v2 p2d_to_netcdf


