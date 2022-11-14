Simulation example of a lateral flux source on a flat topography.
The inlet conditions are given as input of the script. Particles sedimentation and air entrainment are modeled (settling_flag=T, entrainment_flag=T).

A simple firction model, as in [1], is used for this example.
Please note that for this kind of inlet an initial supercritical regime (Ri<1) is required, becuase no characteristic speed should travel upstream.

A Python script (create_example.py) is provided to create the input file for this example (required package: numpy). 
Please provide six arguments:

    print('1) Number of cells\n')
    print('2) Source extent (>0)\n')
    print('3) Initial mfr (>0)\n')
    print('4) Initial temperature (>0)\n')
    print('5) Initial Richardson number (<1)\n')
    print('6) Solid mass fraction (0,1)\n')

1) Number of cells in the x-direction 

2) Source extent (meters)

3) Initial thickness (meters)
 
4) Initial temperature (Kelvin)

5) Initial Richardson number 

6) Initial solid mass fraction (kg/s)

Usage example of the script:

>> python create_example.py 200 250 10000000 900 0.5 0.01

Once the input file (IMEX_SfloW2D.inp) is created, launch the solver:

>> ../../bin/IMEX_SfloW2D

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleLATERAL.bak

The new file can be plotted with Paraview.


