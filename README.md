# Depth-averaged gas-particles model

Shallow water model for multiphase flow (gas+particles) with density of gas temperature-dependent. 

To compile:

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


The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.

It is possible also to use Docker container with the latest version of the model:

> docker pull demichie/sw_var_dens_model_alpine

LINUX/MAC

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v $PWD:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine p2d_to_netcdf

WINDOWS COMMAND LINE

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v %cd%:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine p2d_to_netcdf

WINDOWS POWERSHELL

Create a folder for your simulation with all the input files and then run the container with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine run

If your simulation produced .asc output files, you can post-process those files to have .png files with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine plot_overlay

If your simulation produced .p_2d output files, you can post-process those files to have a netCDF4 file with:

> docker run -v ${PWD}:/home/user_sw/SW_RUNS -i -t demichie/sw_var_dens_model_alpine p2d_to_netcdf


