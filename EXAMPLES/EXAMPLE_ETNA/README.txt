Simulation example of a pyroclastic avalanche with a Voellmy-Salm rheology.

The initial condition is created with a python pre-processing script. The avalanche is computed as a niche obtained by the intersection of the topograpgy and an ellipsoid, defined by two points (x1,y1) (x2,y2) and a semi-axis. The two points in the x-y plane are the coordinates of the lowest and highest point of the initial volume, which also define two semi-axis of the ellipsoid. The length of the third semi-axis is provided by the user. 
The concavity of the niche is searched automatically to match the desired volume. The script generates both the initial thickness of the avalanche and the modified DEM.
The parameters defining the ellipsoid are in the file input_ellipsoid.py. Here an example of its content:

DEM_folder = './DEM/'
DEM_file = 'Etna2014_crop.asc'

x1 = 500300
y1 = 4177600

x2 = 500430
y2 = 4177500

semi_width = 150.0 

vol = 250000

If you need to clean the folder from files created by a previous simulation:

>> ./cleanfolder.sh

To create the initial files, first unzip the file in the DEM folder, and then run:

>> python create_input_ellipsoid.py

To run the example, launch (this assumes that the example is in the original folder):

>> ../../bin/IMEX_SfloW2D

A script to convert the output of the model in NetCDF4 format is provided in the UTILS folder.
First create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/p2d_to_netCDF4.py .

Then, execute the python script

>> python p2d_to_netCDF4.py exampleTAAL.bak

The new file can be plotted with Paraview.

The simulation also creates output files in the ESRI ascii raster format .asc. These files can be used with a GIS software, or post-processed with a python script creating bitmap images.
To do that, first create a simbolic link of the script in this folder:

>> ln -s ../../UTILS/plot_overlay.py .

Then, execute the python script

>> python plot_overlay.py
