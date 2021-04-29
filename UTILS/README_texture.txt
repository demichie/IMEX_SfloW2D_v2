Steps to get a satellite image with the same extent of the simulation domain, to use as texture with Paraview. 
This is based on QGIS, and it requires the plugin QuickMapServices (Plugins > Manage and Install Plugins), and the contributed pack of the plugin (Web > QuickMapServices > Settings > More services > Get contributed pack)
In order to align properly the DEM and the google satellite map, you need to know the UTM zone. You can use this webpage:

https://mangomap.com/robertyoung/maps/69585/what-utm-zone-am-i-in-#


STEPS IN QGIS:
1. Open QGIS
2. Add Raster Layer from panel or Layer > Add Layer > Add Raster Loyer  or Ctrl+Shift+R ('dem_esri.asc')
3. Click on the bottom-right and select the correct Reference Coordinate System (CRS) for the projection, or change it from: Project > Properties. You need to select the appropriate "WGS84 /UTM zone". Once you have selected the UTM zone, click on "Apply" and then "OK".
4. Add the Google Satellite layer: Web > Quick Map Service > Google > Google.cn Satellite
5. Select the Satellite layer
6. Convert map to raster: Processing > Toolbox > Raster tools > Convert map to raster
7. Click the 3 dots for Minimum extent to render and select "Use layer extent", select the raster layer 'dem_esri'
8. Change the resolution. For example, set "Map units per pixel" to 10.0.
9. Click on "Run"
10. Select Output layer
11. Crop to raster extent: Raster > Extraction > Clip Raster by Extent > Use Layer Extent
12. Project > Import/Export > Export Map to Image > (eventually change the dpi ) > Calculate from Layer 'dem_esri' > Save

STEPS IN PARAVIEW:
1. Open Paraview
2. Load the netCDF file (*.nc)
3. Apply the filter "WarpByScalar" to the Scalar "b" to create the 3D topography.
4. Click on the small "2D" icon on the top of the panel with the render view. It should switch to "3D".
5. Select "WarpByScalar1" in the Pipeline Browser, then select the filter "Texture Map to Plane", then click on "Apply".
6. In the Properties tab for "TextureMaptoPlane1", click the wheel to toggle on the advanced properties.
7. Scroll down and, at the beginning of the "Miscellaneous" section, select from the Texture drop-down menu "Load...". Select the texture file saved with QGIS.


