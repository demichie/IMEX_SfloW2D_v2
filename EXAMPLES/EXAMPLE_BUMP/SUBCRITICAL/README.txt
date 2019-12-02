This example simulate a subcritical flow (Ri>1) (see [1], test 3.1.3). 
There is only one phase with constant density (1.D3 kg/s) and no gas. Entrainment, erosion and deposition are neglected in this example. 
No friction is considered in this test (RHEOLOGY_FLAG=F).
On the left boundary (subcritical) only the volumetric flow is fixed: hu = 4.42m^2/s. 
On the right boundary (subcritical) flow height is fixed: h = 2m.

A Python script is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Flow temperature (Kelvin)

Usage example of the script:

>> ./create_example.py 200 1.0 300

Once the input file (SW_VAR_DENS_MODEL.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleSub_200_0024.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleSub_200 100


REFERENCES

[1] Delis, A. I.; Guillard, H. & Tai, Y.-C.
Numerical simulations of hydraulic jumps with the Shear Shallow Water model 
SMAI-Journal of computational mathematics, 2018 , 4 , 319-344
