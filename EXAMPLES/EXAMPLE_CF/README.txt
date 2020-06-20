Simulation example of a collapse of a dilute gas-particle mixture over a topography (Campi Flegrei area). The initial conditions are similar to those presented in [1]

To run the example, unzip first the topography file.

>> unzip topo.zip

Once the topography file is unzipped, create a simbolic link of the executable in this folder:

>> ln -s ../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

Several output files are created as ESRI ascii files (*.asc) and they can be plotted with a GIS.
In addition, model output are saved in other files (*.p2d) which can be converted in netCDF format, and plotted with Paraview. A conversion utility is provided in the UTILS folder, both as a Python script and a Fortran 90 code. They both requires netCDF libraries.

[1] A fast, calibrated model for pyroclastic density currents kinematics and hazard
TE Ongaro, S Orsucci, F Cornolti
Journal of Volcanology and Geothermal Research 327, 257-272
