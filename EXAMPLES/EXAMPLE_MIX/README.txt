This example simulate a 1d flow  with friction, deposition  and entrainment over a gentle slope. 
The total energy conservation equation is solved in this example (ENERGY_FLAG=T).
Topography does not change with deposition (TOPO_CHANGE_FLAG=F), but it is interesting to test the effect on the flow changing the flag.

The settings are similar to those presented in Bursik & Woods 1996 [1]. Here, in addition, the flow has multiple gas components (air and water vapor). The number of additional gas components is defined in the namelist NEWRUN_PARAMETERS (N_ADD_GAS). The parameters for the additional gas components are defined in the namelist GAS_TRANSPORT_PARAMETERS (SP_HEAT_G and SP_GAS_CONST_G). 

A Python script is provided to create the input file for this example. 
Please provide five arguments:

1) Number of cells in the x-direction 

2) Inlet volume fraction of particles

3) Inlet volume fraction of water vapor

4) Inlet flow temperature (Kelvin)

5) Logical for plot of initial solution (true or false)

Usage example of the script:

>> ./create_example.py 400 0.1 0.5 900 false

Run the solver (this assumes that the example is in the original folder):

>> ../../bin/SW_VAR_DENS_MODEL

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py exampleMIX_400_0100.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleMIX_400 100


REFERENCES

[1] Bursik, M. I. &amp; Woods, A. W.
The dynamics and thermodynamics of large ash flows
Bulletin of Volcanology, 1996, 58, 175-193 

