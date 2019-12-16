# Depth-averaged gas-particles model

Shallow water model for multiphase flow (gas+particles) with density of gas temperature-dependent. 

To compile:

> ./configure

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: -fopenmp
2) with intel: -qopenmp

> make

> make install


The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.


