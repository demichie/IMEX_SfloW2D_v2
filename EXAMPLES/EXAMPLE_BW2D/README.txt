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

Usage example of the script:

>> python create_example.py 300 30000.0 2000.0 2.0 330.0 100.0 0.00008 2.5e-4

This is for: domain of 60000mx60000m and 300x300 cells; source radius of 2000m; slope 2 degree; initial mixture temperature 330K; influx velocity 100m/s; initial solid volume fraction 8e-5; particle size 2.5e-4m.

Once the input file (SW_VAR_DENS_MODEL.inp) is created create a simbolic link of the executable in this folder:

>> ln -s ../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py example_BW2D.bak

The new file can be plotted with Paraview.

REFERENCES

[1] Bursik, M. I. & Woods, A. W. The dynamics and thermodynamics of large ash flows. Bulletin of Volcanology, 1996, 58, 175-193 


