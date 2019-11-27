Plot with ParaView 

Load the .nc file
In the Properties box, uncheck "Spherical Coordinates" and click on "Apply".

The following variables are loaded:
T 	temperature (Kelvin)
W	free surface (h+b, in meters)
alphas	volume fraction of solid in the flow
b	topography hight (meters)
dep     deposit thickness (meters)
h	flow thickness (meters)
rhom	flow density (kg/m3)
ux	W-E velocity (m/s)
uy	S-N velocity (m/s)


Apply the filter "WarpByScalar" to the Scalar "b" to create the 3D view of the topography.
Click on the small "2D" icon on the top of the panel with the render view. It should switch to "3D"
On the Color bar, Click on "Solid Color" and change it "b". Click now on the "Edit Color Map" icon, it should open the "Color Map Editor" panel. Choose the colormap you want by clicking on the small icon with a folder and a heart. Close the "Color Map Editor" panel.

In the Pipeline Browser, select the .nc file and apply the filter "Warp by Scalar" to the variable "W" to create the 3D view of the flow over the topography.

In the Pipeline Browser, select the second "WarpByScalar" and apply the "Threshold" filter to remove the areas without the flow. The Scalar "T" can be used for the threshold (increase the Minimum slightly above the atmosphere value, ex. from 300 to 300.01). 
On the Color bar, Click on "Solid Color" and change it "rhom" or "alphas".
Click now on the "Edit Color Map" icon, it should open the "Color Map Editor" panel. Choose the colormap by clicking on the small icon with a folder and a heart. 
Once you have selected the colormap, in the Color Map Editor panel click on "Enable opacity mapping for surfaces". In the Table "Opacity tranfer function values", change the values for "Opacity".
Close the "Color Map Editor" panel.

In the Pipeline Browser select Threshold and apply the "Calculator" filter to create the velocity vector. Change the Result name from "Result" to "Vel". In the box immediately below type:

ux*iHat + uy*jHat

and click on "Apply".

In the Pipeline Browser select Calculator and apply the Glyph filter to it. In the Properties of the Glyph filter, under "Scaling", select "vector" and click on the small reload icon immediately below on the right. Click  on "Apply" on the top of the Properties panel. If no arrows appear, in the "Masking" section of the panel, change from "Uniform Spatial Distribution" to "Every Nth Point". In the "Coloring" section, change from "Solid Color" to "GlyphVector".

In the Pipeline Browser, select the .nc file and apply again the filter "Warp by Scalar" to the variable "b" to create a second the 3D view of the topography.
In the "Coloring" section of the Properties panel, select "dep" 

In the Pipeline Browser, select the last "WarpByScalar" and apply the "Threshold" filter using the variable "dep" to plot only the areas with depost thickness larger than the desired threshold.

From the Filters menu, select AnnotateTimeFilter. In the Properties panel, change the format to "Time: %5.2fs" and click on Apply. 


