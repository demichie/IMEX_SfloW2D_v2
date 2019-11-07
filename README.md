# Welcome to SW_VAR_DENS_MODEL

SW_VAR_DENS_MODEL is a FORTRAN90 code designed to model shallow gas-particles flows over digital elevation models (DEMs) of natural terrain. The model solves for the conservation equations (mass,momentum,energy) of the mixture, and acocunts for sedimentation, erosion, friction and entrainment. The system is described by an hyperbolic system of partial differential equations with relaxation and source terms. It is possible to select a simpler rheology in order to mimic the system of equations described in Kurganov and Petrova, 2007.

Several examples (1D and 2D) are provided with the package.
The code can deal with different scenarios, but its first aim is to treat gravitational flows over topographies described as digital elevation models (DEMs) in the ESRI ascii format. Moreover it is possible to save the solution as as ESRI ascii files, suitable for GIS softwares.

### Authors and Contributors

Mattia de' Michieli Vitturi (@demichie)

### Installation and execution

Check first if you have the LAPACK library installed on your system.

Download the IMEX_SfloW package and create the executable with the following commands from a terminal:

>./configure
>
>make
>
>make install

This will create the executable and copy it in the bin folder. You can test the executable copying it in the EXAMPLES folder and running it.

### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/SW_VAR_DENS_MODEL/wiki](https://github.com/demichie/SW_VAR_DENS_MODEL/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/SW_VAR_DENS_MODEL/html/](http://demichie.github.io/SW_VAR_DENS_MODEL/html/) 

### Acknowledgments

The development of IMEX-SfloW2D has been partially funded by Istituto Nazionale di Geofisica e Vulcanologia.

### References

