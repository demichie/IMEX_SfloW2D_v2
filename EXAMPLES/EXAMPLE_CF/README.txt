Simulation example of a collapse of a dilute gas-particle mixture over a topography (Campi Flegrei area). The initial conditions are similar to those presented in [1]

To run the example, unzip first the topography file.

>> unzip topo.zip

Once the topography file is unzipped, launch the solver:

>> ../../bin/IMEX_SfloW2D

Several output files are created as ESRI ascii files (*.asc) and they can be plotted with a GIS.

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleCF.bak

The new file can be plotted with Paraview.


[1] A fast, calibrated model for pyroclastic density currents kinematics and hazard
TE Ongaro, S Orsucci, F Cornolti
Journal of Volcanology and Geothermal Research 327, 257-272
