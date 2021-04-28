Steps to get a satellite image with the same extent of the simulation domain, to use as texture with Paraview. 
This is based on QGIS, and it requires the plugin QuickMapServices (Plugins > Manage and Install Plugins), and the contributed pack of the plugin (Web > QuickMapServices > Settings > More services > Get contributed pack)

STEPS:
1. Open QGIS
2. Load raster ('dem_esri.asc')
3. Click on the bottom-right and select the correct Reference Coordinate System (CRS) for the projection
4. Add the Google Satellite layer: Web > Quick Map Service > Google > Google.cn Satellite
5. Select the Satellite layer
6. Convert map to raster: Processing > Toolbox > Raster tools > Convert map to raster
7. Click the 3 dots for Minimum extent to render and select "Use layer extent", select the raster layer
8. Change the resolution. For example, set "Map units per pixel" to 10.0.
9. Click on "Run"
10. Crop to raster extent: Raster > Etraction > Clip Raster by Mask Layer
11. Select the raster as "Mask layer" 
12. Project > Import/Export > Export Map to Image
