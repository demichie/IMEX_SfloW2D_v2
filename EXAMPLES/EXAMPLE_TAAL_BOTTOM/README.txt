Simulation example of a flow of a dilute gas-particle mixture over a topography (Taal Volcano complex, Philippines). The equation for pore pressure is solved for this case. 
The radial flow consists of several pulses and then stop. The source is an elliptical area (with semi-axis R and R2 and angle ANGLE_SOURCE), with the gas-particle mixture entering with a vertical velocity from this area. 
The conditions are defined in the following namelist:

&RADIAL_SOURCE_PARAMETERS
 X_SOURCE = 282345.0D0 ,
 Y_SOURCE = 1548428.0D0 ,
 R_SOURCE = 1500.0D0 ,
 R2_SOURCE = 500.0D0 ,
 ANGLE_SOURCE = 135.0 ,
 VEL_SOURCE = 10.0D0 ,
 T_SOURCE = 350.0D0 ,
 ALPHAS_SOURCE = 0.02D0 , 
 TIME_PARAM = 30.D0 , 20.D0 , 0.0D0, 150.D0 ,

The last four parameters define the pulsating flow:
1) the first parameter defines the period between two cycles (30sec) 
2) the second parameter defines the duration of radial flow (20sec). This means that there are 10sec of pause between two pulses.
3) the third parameters defines the time required to the radial source to increase from no flow to maximum flow, and vice versa.
4) the last parameter is the total duration of radial source.

For a steady state source, remove the line with the four parameters.

For a single pulse, the first parameter (time between cycles) has to be larger than the fourth one (duration of source). For example:

TIME_PARAM = 300.D0 , 20.D0 , 0.0D0, 150.D0 ,

In this case, we have a single 20 seconds pulse.


To run the example, first create a simbolic link of the executable in this folder:

>> ln -s ../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py

The new file can be plotted with Paraview.
