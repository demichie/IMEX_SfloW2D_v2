This example simulate a transcritical flow (from Ri>1, subcritical, to Ri<1, supercritical) without shocks (see [1], test 3.1.4). 
There is only one phase with constant density (1.D3 kg/s) and no gas. Entrainment, erosion and deposition are neglected in this example. 
No friction is considered in this test (RHEOLOGY_FLAG=F).
On the left boundary (subcritical) only the volumetric flow is fixed. 
On the right boundary (supercritical, right-going) zero gradient is fixed for all variables.

A Python script is provided to create the input file for this example. 
Please provide six arguments:

1) Number of cells

2) Variables to reconstruct: phys or cons

3) Order of the RK scheme

4) Initial solid volume fraction (0-1)

5) Erosion coefficient (>0)

6) Settling velocity (>0)

7) Temperature

Usage example of the script:

>> ./create_example.py 200 phys 2 1.0 0 0 300

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleTrShock_200_0024.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleTrShock_200 100


REFERENCES

[1] Delis, A. I.; Guillard, H. & Tai, Y.-C.
Numerical simulations of hydraulic jumps with the Shear Shallow Water model 
SMAI-Journal of computational mathematics, 2018 , 4 , 319-344
