Example of ensemble simulation for pyroclastic avalanches at Etna. 
For ensamble scenarios, in order to avoid multiple copies of the same files, we suggest to keep the structure of the folders as in this example. 

-- ENSEMBLE FOLDER
      |
      |---- DEM
      |
      |---- templatedir
                 |
                 |---- DEM
                 
The executable should be placed in the templatedir as a symbolic link.                 
The topography file should be in the top DEM folder, and a symbolic link should be created in the templatedir/DEM folder.

STEP 1: create a symbolic link to the executable in the template folder

> cd templatedir
> ln -s ../../../bin/IMEX_SfloW2D .
> cd ..

STEP 2: unzip the DEM file and create a symbolic link in the templatedir/DEM folder 

> cd DEM
> unzip Etna2014_crop.zip
> cd ../templatedir
> mkdir DEM
> cd DEM
> ln -s ../../DEM/Etna2014_crop.asc .
> cd ../..

STEP 3: generate the input parameters of the ensemble, creating a CSV file

> python create_ensemble.py

The number of sample, the uncertain parameters and their ranges are defined in the file create_ensemble.py.
Each parameters has a keyword, for example "var1".
Here, we have two different locations with different weights and uncertainty ranges for the rheological parameters.

STEP 4: create the folders for the ensemble

> python create_inputfiles.py

This script create the folders by making a copy of the folder "templatedir". The values of the uncertain parameters
in the input files must be replaced by the strings "ENSAMBLE_xyz", where "xyz" is the keyword assigned in create_ensamble.py. For example, if "var1" is defined as unknown in create_ensemble.py, we must have in an input file the string "ENSEMBLE_var1".   

In the is simulation we have, for the file IMEX_SfloW2D.inp:

&RESTART_PARAMETERS
 RESTART_FILE ="./DEM/init_volume.asc" ,
 T_INIT = 373.D0 ,
 T_AMBIENT = 300.D0 ,
 SED_VOL_PERC = 50.0  ,
 u_init = ENSEMBLE_vx,
 v_init = ENSEMBLE_vy,
/ 

&RHEOLOGY_PARAMETERS
 RHEOLOGY_MODEL = 1,
 MU = ENSEMBLE_mu,
 XI = ENSEMBLE_xi
 /

STEP 5: launch the simulations

> python launch_jobs.py

This script launch all the simulations. The maximum number of concurrent runs can be prescribed in the file (max_processes). This script launches the batch script "run.sh" present in each folder "ensemble.xxxxx". This script executes the pre-processing step for the initial volume, run the code and post-processes the output files. If the post-processing is not needed, comment the lines in "run.sh" inside tmeplatedir (before creating the folders).

STEP 6: create probabilistic output

> python create_prob_maps.py

This script creates probabilistic maps (as raster .asc files) from the .asc files for the maximum thickness and dyn.press.+thickness (VT files). The elements of the ensamble can be treated as a single large group or partitioned into smaller groups by using the pandas "groupby" function. A few examples can be found in the script. For each group, a folder named result.xxx is created.

STEP 8: post-processing of probabilistic output

> python plot_ensemble.py

This script loops over the different folders result.xxx, and creates .png maps for each probabilistic output. 


 

