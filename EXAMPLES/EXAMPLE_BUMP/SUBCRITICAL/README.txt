This example simulate a subcritical flow (Ri>1) (see [1], test 3.1.3). 
There is only one phase with constant density (1.D3 kg/s) and no gas. Entrainment, erosion and deposition are neglected in this example. 
No friction is considered in this test (RHEOLOGY_FLAG=F).
On the left boundary (subcritical) only the volumetric flow is fixed: hu = 4.42m^2/s. 
On the right boundary (subcritical) flow height is fixed: h = 2m.

A Python script is provided to create the input file for this example. 
Please provide two arguments:

1) Number of cells in the x-direction 

2) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 200 false

Run the solver (this assumes that the example is in the original folder):

>> ../../../bin/IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> python plot_phys.py exampleSub_200_0024.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> python plot_animated.py exampleSub_200 100


REFERENCES

[1] Delestre, O., Lucas, C., Ksinant, P. A., Darboux, F., Laguerre, C., Vo, T. N. T., ... & Cordier, S. (2013). SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. International Journal for Numerical Methods in Fluids, 72(3), 269-300.
