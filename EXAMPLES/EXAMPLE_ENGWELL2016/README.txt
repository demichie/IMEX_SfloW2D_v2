Simulation example of a radial flux source on a flat topography. This example is aimed at reproducing the results of Engwell at al. [2], Figure 3.
The inlet conditions are given as input of the script. Particles sedimentation and air entrainment are modeled (settling_flag=T, entrainment_flag=T).

A simple firction model, as in [1], is used for this example.

A Python script (create_example.py) is provided to create the input file for this example (required package: numpy). 
Please provide six arguments:

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

5) Initial Richardsson number

6) Initial solid mass fraction

Usage example of the script:

>> python create_example.py 200 2000 2000 900 0.9 0.8

Once the input file (SW_VAR_DENS_MODEL.inp) is created, create a simbolic link of the executable in this folder:

>> ln -s ../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleRS.bak

The new file can be plotted with Paraview.

REFERENCES

[1] Bursik, M. I. & Woods, A. W. The dynamics and thermodynamics of large ash flows. Bulletin of Volcanology, 1996, 58, 175-193 
[2] Engwell, S.L., de' Michieli Vitturi, M. Esposti Ongaro, T,. and Neri, A. Insights into the formation and dynamics of coignimbrite plumes from one-dimensional models. Journal of Geophysical Research: Solid Earth 2016, 121 (6), 4211–4231


