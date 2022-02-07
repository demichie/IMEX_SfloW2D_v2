Simulation example of a radial inlet flux source on a flat topography. The inlet conditions are given as input of the script. Particles sedimentation and air entrainment are not considered (settling_vel=0, entrainment_flag=F).

A simple firction model, as in [1], is used for this example.

A Python script (create_example.py) is provided to create the input file for this example (required packages: numpy, basemap). 
Please provide three arguments:

    print('1) Number of cells\n')
    print('2) Source radius (>0)\n')
    print('3) Initial thickness (>0)\n')
    print('4) Temperature (>0)\n')
    print('5) Radial velocity (>0)\n')
    print('6) Solid volume fraction (0,1)\n')

1) Number of cells in the x-direction 

2) Source radius

3) Initial thickness

4) Initial temperature

5) Initial radial velocity

6) Initial solid volume fraction

Usage example of the script:

>> python create_example.py 150 20 0.08 300 2.5 0.99

Once the input file (SW_VAR_DENS_MODEL.inp) is created, run the solver (this assumes that the example is in the original folder):

>> ../../bin/SW_VAR_DENS_MODEL

Two Python scripts to plot the results are provided. The first one works better with small simulations.
Usage example:

>> ./plot_small.py exampleRS_0050.p_2d

The second script allows to plot the topography, the flow thickness and the values of an additional variable (required packages: pandas, plotly). The plot is shown in a browser.

Usage example:

>> ./plot_large.py exampleRS_0050.p_2d u

In addition, a script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleRS.bak

The new file can be plotted with Paraview.


REFERENCES

[1] Bursik, M. I. & Woods, A. W. The dynamics and thermodynamics of large ash flows. Bulletin of Volcanology, 1996, 58, 175-193 


