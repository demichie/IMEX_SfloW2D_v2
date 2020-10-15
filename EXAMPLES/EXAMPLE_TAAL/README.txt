Simulation example of a supercritical radial flow of a dilute gas-particle mixture over a topography (Taal Volcano complex, Philippines). The radial flow consists of several pulses and then stop. The conditions are defined in the following namelist:

&RADIAL_SOURCE_PARAMETERS
 X_SOURCE = 282345.0D0 ,
 Y_SOURCE = 1548428.0D0 ,
 H_SOURCE = 40.0D0 ,
 R_SOURCE = 250.0D0 ,
 VEL_SOURCE = 50.0D0 ,
 T_SOURCE = 350.0D0 ,
 ALPHAS_SOURCE = 0.02D0 , 
 TIME_PARAM = 30.D0 , 20.D0 , 4.0D0, 150.D0 ,

The last four parameters define the pulsating flow:
1) the first parameter defines the period between two cycles (30sec) 
2) the second parameter defines the duration of radial flow (20sec). This means that thera are 10sec of pause between two pulses.
3) the third parameters defines the time required to the radial source to increase from no flow to maximum flow, and vice versa. Initial thickness and velocity are scaled in order to keep constant the Richardson number.
4) the toal duration of radial source.

To run the example, first create a simbolic link of the executable in this folder:

>> ln -s ../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleTAAL.bak

The new file can be plotted with Paraview.
