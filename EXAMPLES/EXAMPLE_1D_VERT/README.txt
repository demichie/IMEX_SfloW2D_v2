This example simulate a 1D flow  with friction, deposition  and entrainment over a gentle slope. The temperature equation is solved in this example (ENERGY_FLAG=F).

The vertical profiles for particle concentration and velocity can be activated by setting, in the NEWRUN_PARAMETERS namelist:

>> VERTICAL_PROFILES_FLAG = T

Please note that the mass flow rates at the inlet, with and without the vertical profiles, are different.

A Python script is provided to create the input file for this example. 
Please provide four arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Flow temperature (Kelvin)

4) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 400 0.01 900 false

Run the solver (this assumes that the example is in the original folder):

>> ../../bin/IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleBW_400_0100.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleBW_400 100

