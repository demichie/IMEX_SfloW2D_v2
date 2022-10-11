This example simulates a 2d radial flow  with friction, deposition and entrainment over a slope. The settings are similar to those presented in Bursik & Woods 1996 [1].
The total energy conservation equation is solved in this example (ENERGY_FLAG=T).
Topography does not change with deposition (TOPO_CHANGE_FLAG=F). A simple firction model, as in [1], is used for this example.

Here, instead of a radial source, flow enters from the bottom of a circular area centered at (0,0), without horizontal momentum. Deposition and entrainment are not present in this area. In this way, when a steady condition is reached, the mass flux entering the domain from the bottom source is the same of the mass flux going out of the source area. 

A Python script (create_example.py) is provided to create the input file for this example (required packages: numpy, basemap). 
Please provide six arguments:

1) Number of cells in the x-direction 

2) Half-domain size (m)

3) Source radius (m)

4) Topography slope (degree)

5) Initial temperature (K)

6) Inflow velocity (m/s)

7) Initial solid volume fraction

8) Particle size (m)



The last four parameters define the pulsating flow:
1) the first parameter defines the period between two cycles (30sec) 
2) the second parameter defines the duration of radial flow (20sec). This means that thera are 10sec of pause between two pulses.
3) the third parameters defines the time required to the radial source to increase from no flow to maximum flow, and vice versa. Initial thickness and velocity are scaled in order to keep constant the Richardson number.
4) the toal duration of radial source.


Usage example of the script:

>> python create_example.py 300 30000.0 2000.0 2.0 330.0 100.0 0.00008 2.5e-4

This is for: domain of 60000mx60000m and 300x300 cells; source radius of 2000m; slope 2 degree; initial mixture temperature 330K; influx velocity 100m/s; initial solid volume fraction 8e-5; particle size 2.5e-4m.

With these values, the following namelist is created in the input file:

&RADIAL_SOURCE_PARAMETERS
 X_SOURCE=  0.0D0     ,
 Y_SOURCE=  0.0D0     ,
 R_SOURCE= 2000.0D0     ,
 VEL_SOURCE= 100.0D0     ,
 T_SOURCE= 330.0D0     ,
 ALPHAS_SOURCE= 8e-05       ,
 TIME_PARAM = 600.D0 , 100.D0 , 50.0D0, 50.D0 ,
 /

The last four parameters define the temporal changes in the inflow:
1) the first parameter defines the period between two cycles (600sec)
2) the second parameter defines the duration of a source pulse (100sec). This means that there are 500sec of pause between two pulses.
3) the third parameters defines the time required to the radial source to change linearly from no flow to maximum flow, and vice versa. In this case, this time is half of the pulse duration, i.e. source increase linearly for 50s and then immediately decrease linearly.
4) the point within a cycle corresponding to the beginning of the simulation. In this case we start at 50s of a cycle, i.e. at the time corresponding to the peak. With these four values of the parameters, the source starts from the peak, decreases linearly for 50s, and then remains null for the remaining of the simulation.


Finally, launch the solver:

>> ../../bin/IMEX_SfloW2D

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py example_BW2D.bak

The new file can be plotted with Paraview.

REFERENCES

[1] Bursik, M. I. & Woods, A. W. The dynamics and thermodynamics of large ash flows. Bulletin of Volcanology, 1996, 58, 175-193 


