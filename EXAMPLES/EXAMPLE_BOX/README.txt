This example simulate a 1D dam-break problem with erosion and deposition [1] over a flat topography. Erosion and deposition are proportional to two input parameters: an empirical erosion coefficient (EROSION_COEFF) and the sediment settling velocity (SETTLING_VEL). 
No friction is considered in this test (RHEOLOGY_FLAG=F).

A Python script is provided to create the input file for this example. 
Please provide six arguments:

1) Number of cells

2) Variables to reconstruct: phys or cons

3) Order of the RK scheme

4) Initial solid volume fraction (0;1)

5) Erosion coefficient (>0)

6) Settling velocity (>0)

7) Temperature

Usage example of the script:

>> ./create_example.py 400 phys 2 0.5 0.0 1.0 300

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleErosionDeposition_400_0010.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleErosionDeposition_400 100


REFERENCES

[1] Fagents, S. A., & Baloga, S. M. (2006). Toward a model for the bulking and debulking of lahars. Journal of Geophysical Research: Solid Earth, 111(B10).
