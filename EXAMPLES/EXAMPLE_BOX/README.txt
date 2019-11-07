This example simulate a 1D dam-break problem with deposition, entrainment and friction over a flat topography. 

A Python script is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Erosion coefficient ( 0 => no erosion )

4) Sedimentation coefficient ( 0 => no deposition )

5) Flow temperature (Kelvin)

Usage example of the script:

>> ./create_example.py 400 0.5 0.0 1.0 300

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleBOX_400_0010.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleBOX_400 100
<<<<<<< HEAD


REFERENCES

