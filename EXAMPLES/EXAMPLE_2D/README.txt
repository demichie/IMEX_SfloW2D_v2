Simulation example of an avalanche of finite granular mass sliding down an inclined plane and merging continuously into a horizontal plane is presented. The initial conditions and the topography of this tutorial are similar to those in Example 4.1 from Wang et al., 2004 [1]. A hemispherical shell holding the material together is suddenly released so that the bulk material commences to slide on an inclined flat plane at 35° into a horizontal run-out plane connected by a smooth transition.

Please note that we do not use the same rheological model as in the original paper of the example.

A Python script (create_example4.py) is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells in the x-direction 

2) Volume fraction of particles

3) Erosion coefficient ( 0 => no erosion )

4) Sedimentation coefficient ( 0 => no deposition )

5) Flow temperature (Kelvin)

Usage example of the script:

>> python create_example.py 100 0.5 0.0 1.0 300

Once the input file (SW_VAR_DENS_MODEL.inp) is created create a simbolic link of the executable in this folder:

>> ln -s ../../bin/SW_VAR_DENS_MODEL .

Finally, launch the solver:

>> ./SW_VAR_DENS_MODEL

A Python script to plot the results is provided. With this script you can choose the output to plot.
Usage example:

>> ./plot_p_2d.py example2D_0030.p_2d

A Python script to create an animation of the simulation is also provided. This script plot the topography and animate the flow over it. An animated gif is saved by the script. The script requires two argument: run name (as given in SW_VAR_DENS_MODEL.inp) and duration of the animated gif in seconds

Usage example:

>> ./animate_asc.py example2D 10


REFERENCES

[1] Wang, Y., K. Hutter and S. P. Pudasaini, The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud, ZAMM · Z. Angew. Math. Mech. 84, No. 8, 507 – 527 / DOI 10.1002/zamm.200310123, 2004.

