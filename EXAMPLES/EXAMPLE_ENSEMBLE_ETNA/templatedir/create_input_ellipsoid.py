#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%
% Pre-processing tool for IMEX_SfloW2D simulations
%
% by M. de' Michieli Vitturi (2021-10-14)
%
% This script generate an .asc file with the initial thickness for the simulation,
% obtaining by subtracting from the original den an ellipsoidal niche.
% The following input parameters are required:
% - DEM files (ESRII ascii format)
% - x,y, coordinates of highest point P1(UTM coordinates)
% - x,y, coordinates of highest point P2 (UTM coordinates)
% - horizontal semiaxis of the ellipsoid orthogonal to the direction from P1 to P2 
% The following output files are saved:
% - modified DEM
% - initial landslide volume
% 
"""
import numpy as np                      
from linecache import getline
import sys
import os.path
import time
import re
import pandas as pd

print('')
print('----------------------------------------------')
print("Pre-processing tool by M. de' Michieli Vitturi")
print('----------------------------------------------')
print('')


# Parameters are read from file input.py
from input_ellipsoid import *

# Function to convert   
def listToString(s):  
    
    # initialize an empty string 
    str1 = ""  
    
    # traverse in the string   
    for ele in s:  
        str1 += str(ele)+" "   
    
    # return string   
    return str1  

def isfloat(string):
  try:
    return float(string) and '.' in string  # True if string is a number contains a dot
  except ValueError:  # String is not a number
    return False

def compute_vol_new(DEM,dist,XS,YS,x1,y1,semi_axis,aspect_ratio,cell_topo,n):

    indexes = np.where( ( dist>semi_axis-cell_topo )*( dist<semi_axis+cell_topo ) )
    indexes_np = np.asarray(indexes)

    """
    idx_max = np.argmax(DEM[indexes])

    ix = indexes_np[1,idx_max]
    iy = indexes_np[0,idx_max]
    """
    
    idx_min = np.argmin(DEM[indexes])

    ix = indexes_np[1,idx_min]
    iy = indexes_np[0,idx_min]

    dist_real = np.sqrt( (XS[iy,ix]-x1)**2 + (YS[iy,ix]-y1)**2 )

    x2 = x1 + (XS[iy,ix]-x1)*semi_axis/dist_real
    y2 = y1 + (YS[iy,ix]-y1)*semi_axis/dist_real

    # compute the first point elevation by linearly interpolating the DEM data
    X1 = np.zeros((1,1))
    X1[0,0] = x1

    Y1 = np.zeros((1,1))
    Y1[0,0] = y1

    z1 = interp2Dgrids(xs_topo,ys_topo,DEM,X1,Y1).flatten()[0]

    # compute the second point elevation by linearly interpolating the DEM data
    X2 = np.zeros((1,1))
    X2[0,0] = x2

    Y2 = np.zeros((1,1))
    Y2[0,0] = y2

    z2 = interp2Dgrids(xs_topo,ys_topo,DEM,X2,Y2).flatten()[0]

    # define the center of the ellipsoid 3D coordinates and the vertical semi-axis
    if z1>z2:

        xc = x2
        yc = y2
        zc = z1
    
        xh = x1
        yh = y1
        c = z1-z2
    
    else:

        xc = x1
        yc = y1
        zc = z2
    
        xh = x2
        yh = y2
        c = z2-z1
    
    print('Ellissoid center:',xc,yc,zc)    

    semi_width = semi_axis

    # define the other semi-axis of the ellipsoid
    a = np.sqrt((y1-y2)**2+(x1-x2)**2)
    b = semi_width / aspect_ratio

    print('Ellipsoid semi-axis: ',a,b,c)

    # compute the horizontal angle between the original x,y, coordinate system
    # and the coordinate system defined by the horizonal semi-axis a,b
    alpha_deg = -np.arctan2(yh-yc,xh-xc)/np.pi*180.0

    print('Rotation angle: ',alpha_deg)

    # translate the coordinate of the centers of the cells
    X,Y = np.meshgrid(xs_topo-xc,ys_topo-yc)

    # rotation to find the coordinates of the centers of the cells
    # in the new coordinate system, centered in the center of the
    # ellispoid and with the axis parallel to its semi-axis
    alpha_rad = alpha_deg/180.0*np.pi

    X_new = X*np.cos(alpha_rad) - Y*np.sin(alpha_rad)
    Y_new = X*np.sin(alpha_rad) + Y*np.cos(alpha_rad)

    Z_new,Z_ell,vol = compute_vol(X_new,Y_new,DEM,zc,a,b,c,n)
    
    return Z_new,Z_ell,vol,zc,alpha_rad


def compute_vol(X_new,Y_new,DEM,zc,a,b,c,n):

    # compute the ellipsoid
    Z_temp1 = zc - c * ( np.maximum(np.zeros_like(X_new),1.0 - (np.abs(X_new/a))**2 - np.abs((Y_new/b))**2) )**(1.0/n)
    Z_temp2 = Z_temp1 * ( (X_new/a)**2 + (Y_new/b)**2 <1 ) + 1.e5 ** ( (X_new/a)**2 + (Y_new/b)**2 >=1 )

    # compute the modified DEM
    Z_new = np.minimum(DEM,Z_temp2)
    vol = np.sum(DEM-Z_new)*cell_topo**2

    print('Initial volume',np.sum(DEM-Z_new)*cell_topo**2,' m3')

    return Z_new,Z_temp1,vol


def regrid(xin, yin, fin, xl, xr , yl, yr):

        
    nXin = xin.shape[0]-1
    nYin = yin.shape[0]-1

    dXin = xin[1] - xin[0]
    dYin = yin[1] - yin[0]
    
    ix1 = np.maximum( 0 , np.ceil(( xl - xin[0] ) / dXin ).astype(int) -1 )
    ix2 = np.minimum( nXin , np.ceil( ( xr -xin[0] ) / dXin ).astype(int) )
        
    iy1 = np.maximum( 0 , np.ceil( ( yl - yin[0] ) / dYin ).astype(int) -1 )
    iy2 = np.minimum( nYin , np.ceil( ( yr - yin[0] ) / dYin ).astype(int) )

    fout = 0.0

    for ix in range(ix1,ix2):
        
       alfa_x = ( np.minimum(xr,xin[ix+1]) - np.maximum(xl,xin[ix]) ) / ( xr - xl )

       for iy in range(iy1,iy2):
                
          alfa_y = ( np.minimum(yr,yin[iy+1]) - np.maximum(yl,yin[iy]) ) / ( yr - yl )

          fout = fout + alfa_x * alfa_y * fin[iy,ix]

    return fout
 

# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def interp2Dgrids(xin,yin,Zin,Xout,Yout):
    """
    Interpolation from a regular grid to a second regular grid 

    @params:
        xin      - Required : original grid X values (1D Dble) 
        yin      - Required : original grid Y values (1D Dble) 
        Zin      - Required : original grid Z values (2D Dble) 
        xout     - Required : new grid X values (2D Dble) 
        yout     - Required : new grid Y values (2D Dble) 
    """
    xinMin = np.min(xin)
    xinMax = np.max(xin)

    yinMin = np.min(yin)
    yinMax = np.max(yin)

    cellin = xin[1]-xin[0]

    if ( Xout.ndim == 2 ):

        xout = Xout[0,:]
        xoutMin = np.min(xout)
        xoutMax = np.max(xout)
        
    else:     

        xout = Xout
        
    xoutMin = np.min(Xout)
    xoutMax = np.max(Xout)

    if ( Yout.ndim == 2 ):

        yout = Yout[:,0]
        
    else:
    
        yout = Yout
    
    youtMin = np.min(yout)
    youtMax = np.max(yout)


    # Search for the cell containing the center of the parent lobe
    xi = (xout - xinMin)/cellin
    yi = (yout - yinMin)/cellin

    # Indexes of the lower-left corner of the cell containing the center of the parent lobe     
    ix = np.maximum(0,np.minimum(xin.shape[0]-2,np.floor(xi).astype(int)))
    iy = np.maximum(0,np.minimum(yin.shape[0]-2,np.floor(yi).astype(int)))

    # Indexes of the top-right corner of the cell containing the center of the parent lobe     
    ix1 = ix+1
    iy1 = iy+1

    # Relative coordinates of the center of the parent lobe in the cell
    
    xi_fract = np.maximum(0.0,np.minimum(1.0,(xi-ix).reshape(1,Xout.shape[1])))
        
    yi_fract = np.maximum(0.0,np.minimum(1.0,(yi-iy).reshape(Yout.shape[0],1)))
        
    xi_out_yi = np.outer(yi_fract,xi_fract)

    Zout = xi_out_yi * Zin[np.ix_(iy1,ix1)] + \
           ( xi_fract - xi_out_yi ) * Zin[np.ix_(iy,ix1)] + \
           ( yi_fract - xi_out_yi ) * Zin[np.ix_(iy1,ix)] + \
           ( 1.0 - xi_fract - yi_fract + xi_out_yi ) * Zin[np.ix_(iy,ix)] 

    return Zout



# ------------------- READ TOPOGRAPHY ------------------------------------------
# Parse the topography header 
DEM_file = DEM_folder+DEM_file

# Load the topography into a numpy array
print('Reading DEM file: '+DEM_file)

hdr = [getline(DEM_file, i) for i in range(1,7)]
values = [float(h.split(" ")[-1].strip()) \
 for h in hdr]
cols_topo,rows_topo,lx_topo,ly_topo,cell_topo,nd = values
cols_topo = int(cols_topo)
rows_topo = int(rows_topo)

# grid with coordinates at cell centers
xs_topo = lx_topo+0.5*cell_topo + np.linspace(0,(cols_topo-1)*cell_topo,int(cols_topo))
ys_topo = ly_topo+0.5*cell_topo + np.linspace(0,(rows_topo-1)*cell_topo,int(rows_topo))

XS,YS = np.meshgrid(xs_topo,ys_topo)

dist = np.sqrt( (XS-x1)**2 + (YS-y1)**2 ) 


#print('indexes',indexes)
#print(indexes_np.shape)

DEM = pd.read_table(DEM_file, delim_whitespace=True, header=None,skiprows=6, dtype='unicode').astype(float).values
DEM = np.flipud(DEM)

print('Completed reading DEM, shape: '+str(DEM.shape))
print('')



semi_axis0 = 0.5*( 3.0/4.0*np.pi * vol)**(1.0/3.0)
print('semi_axis0',semi_axis0)
        
Z_new_temp,Z_ell_temp,vol0,zc,alpha_rad = compute_vol_new(DEM,dist,XS,YS,x1,y1,semi_axis0,aspect_ratio,cell_topo,2.0)

if vol0 > vol:

    Z_new = Z_new_temp
    Z_ell = Z_ell_temp

else:

    semi_axis2 = 4.0*( 3.0/4.0*np.pi * vol)**(1.0/3.0)
    print('semi_axis2',semi_axis2)
        
    Z_new_temp,Z_ell_temp,vol2,zc,alpha_rad = compute_vol_new(DEM,dist,XS,YS,x1,y1,semi_axis2,aspect_ratio,cell_topo,2.0)


    
if vol2 < vol:    

    Z_new = Z_new_temp
    Z_ell = Z_ell_temp

elif vol0<vol:

    for i in range(20):
            
        semi_axis1 = semi_axis0+((vol-vol0)/(vol2-vol0))**(1.0/3.0)*(semi_axis2-semi_axis0)
        print('semi_axis1',semi_axis1)
        
        Z_new,Z_ell,vol1,zc,alpha_rad = compute_vol_new(DEM,dist,XS,YS,x1,y1,semi_axis1,aspect_ratio,cell_topo,2.0)
        
        if ( vol1 < vol ):
        
            semi_axis0 = semi_axis1
            
        else:

            semi_axis2 = semi_axis1
                        

# Write output
print('')
print('Saving output')

# modified DEM
header = "ncols     %s\n" % Z_new.shape[1]
header += "nrows    %s\n" % Z_new.shape[0]
header += "xllcorner " + str(lx_topo) +"\n"
header += "yllcorner " + str(ly_topo) +"\n"
header += "cellsize " + str(cell_topo) +"\n"
header += "NODATA_value " + str(nd)

np.savetxt(DEM_folder+'new_dem.asc', np.flipud(Z_new), header=header, fmt='%1.5f',comments='')

# intial volume
header = "ncols     %s\n" % Z_new.shape[1]
header += "nrows    %s\n" % Z_new.shape[0]
header += "xllcorner " + str(lx_topo) +"\n"
header += "yllcorner " + str(ly_topo) +"\n"
header += "cellsize " + str(cell_topo) +"\n"
header += "NODATA_value " + str(int(nd))

np.savetxt(DEM_folder+'init_volume.asc', np.flipud(DEM-Z_new), header=header, fmt='%1.5f',comments='')


# complete ellipsoidal
header = "ncols     %s\n" % Z_new.shape[1]
header += "nrows    %s\n" % Z_new.shape[0]
header += "xllcorner " + str(lx_topo) +"\n"
header += "yllcorner " + str(ly_topo) +"\n"
header += "cellsize " + str(cell_topo) +"\n"
header += "NODATA_value " + "{:1.5f}".format(zc)

np.savetxt(DEM_folder+'ellipsoid.asc', np.flipud(Z_ell), header=header, fmt='%1.5f',comments='')

vx = -mod_vel*np.cos(alpha_rad)
vy = mod_vel*np.sin(alpha_rad)

print('vx,vy',vx,vy)
  
# Opening our text file in read only
# mode using the open() function
with open(r'IMEX_sfloW2D.inp', 'r') as file:
  
    # Reading the content of the file
    # using the read() function and storing
    # them in a new variable
    data = file.read()
  
    # Searching and replacing the text
    # using the replace() function
    data = data.replace('ENSEMBLE_vx', str(vx))
    data = data.replace('ENSEMBLE_vy', str(vy))
  
# Opening our text file in write only
# mode to write the replaced content
with open(r'IMEX_sfloW2D.inp', 'w') as file:
  
    # Writing the replaced data in our
    # text file
    file.write(data)
  
# Printing Text replaced
print("Text replaced")
