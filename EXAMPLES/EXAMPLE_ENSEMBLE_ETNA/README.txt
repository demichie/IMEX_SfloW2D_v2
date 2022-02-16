Example of ensemble simulation for pyroclastic avalanches at Etna

STEP 1: create a symbolic link in the template folder to the executable

> cd templatedir
> ln -s ../../../bin/SW_VAR_DENS_MODEL .
> cd ..

STEP 2: unzip the DEM file and create a symbolic link in the templatedir/DEM folder to it

> cd DEM
> unzip Etna2014_crop.zip
> cd ../templatedir/DEM
> ln -s ../../DEM/Etna2014_crop.asc .
> cd ../..

STEP 3: create a CSV file with uncertain input parameters

> python create_ensemble.py

The number of sample, the uncertain parameters and their ranges are defined in the create_ensemble.py file
Each parameters has a keyword, for example "var1"

STEP 4: create the ensamble folders

> python create_inputfiles.py

This script create the folders by making a copy of the folder "templatedir". The values of the uncertain parameters
in the input file must be replaced by the string "ENSAMBLE_xyz", where "xyz" is the keyword assigned in create_ensamble.py. For example, if "var1" is defined as unknown in create_ensemble.py, we must have in an input file "ENSEMBLE_var1".   

STEP 5: launch the simulations

> python launch_jobs.py

This script launch all the simulations. The maximum number of concurrent runs can be prescribed in the file. This script launch a batch script "run.sh" present in each folder "ensemble.xxxxx". If the post-processing is not needed, comment the lines in "run.sh" inside tmeplatedir (before creating the folders).

STEP 6: create probabilistic output

> python create_prob_maps.py

This script creates probabilistic maps (as raster .asc files) from the .asc files for the maximum thickness and dyn.press.+thickness (VT files). The elements of the ensamble can be treated as a single large group or partitioned in smaller groups by using the pandas "groupby" function. A few examples can be found in the script. For each group, a folder named result.xxx is created.

STEP 8: post-processing of probabilistic output

> python plot_ensemble.py

This script loops over the different folders result.xxx, and creates .png maps for each probabilistic output. 


 

