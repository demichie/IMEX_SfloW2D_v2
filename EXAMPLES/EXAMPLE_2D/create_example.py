#!/usr/bin/env python
"""
% This function 

"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import matplotlib.tri as mtri

import time
import sys

if len(sys.argv)==4: 

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        nx_cells = int(a)
 
    else:
 
        sys.exit()

    a = sys.argv[2]
    try:
        alfas = float(a)
    except ValueError:
        print("You must enter a float")

    a = sys.argv[3]
    try:
        T = float(a)
    except ValueError:
        print("You must enter a float for temperature: "+a)

else:

    print('Please provide three arguments:\n')
    print('1) Number of cells\n')
    print('2) Solid volume fraction (0,1)\n')
    print('3) Temperature (>0)\n')
    sys.exit()


n_solid = 1

rho_s = 2000.0
rho_a = 1.2
SP_HEAT_S = 1617.0
SP_HEAT_A = 998.0
SP_GAS_CONST_A = 287.051
PRES = 101300.0

rho_a = PRES / ( T * SP_GAS_CONST_A )

rho_m = alfas * rho_s + ( 1.0-alfas ) * rho_a

xs = alfas * rho_s / rho_m

SP_HEAT_MIX = xs * SP_HEAT_S + ( 1.0 - xs ) * SP_HEAT_A



# Define the boundaries x_left and x_right of the spatial domain
x_min = 0.0
x_max = 30.0

y_min = -10.0
y_max = 10.0


# Define the number n_points of points of the grid
nx_points  = nx_cells+1

# Define the grid stepsize dx
dx = ( x_max - x_min ) / ( nx_cells )

# print('dx',dx,x_max - x_min, nx_cells)

# Define the array x of the grid points
x = np.linspace(x_min,x_max,nx_points)

x_cent = np.linspace(x_min+0.5*dx,x_max-0.5*dx,nx_cells)



dy = dx

# print('dy',dy)
ny_half_cells = int(np.ceil(y_max/dy))
ny_cells = 2*ny_half_cells
ny_points = ny_cells+1

y_min = -dy*ny_half_cells
y_max = -y_min

print('Number of cells in the y-direction')
print(ny_cells)

n_cells = nx_cells * ny_cells

# print(y_min)
# print(y_max) 

# Define the array x of the grid points
y = np.linspace(y_min,y_max,ny_points)

y_cent = np.linspace(y_min+0.5*dy,y_max-0.5*dy,ny_cells)

X, Y = np.meshgrid(x, y)
X_cent, Y_cent = np.meshgrid(x_cent, y_cent)

# print X.shape
# print X_cent.shape

Z = np.zeros_like(X)

Z_cent = np.zeros_like(X_cent)
W_cent = np.zeros_like(X_cent)
H_cent = np.zeros_like(X_cent)
U_cent = np.zeros_like(X_cent)
V_cent = np.zeros_like(X_cent)


# define the topography
for i in range(nx_points-1,-1,-1):

    if X[0,i]>21.5:
    
        Z[:,i] = 0
    
    elif X[0,i]>17.5:
    
        Z[:,i] = Z[:,i+1] + np.sin(np.deg2rad(35.0))*(1.0-0.25*(X[:,i]-17.5))*dx
 
    else:
    
        Z[:,i] = Z[:,i+1] + np.sin(np.deg2rad(35.0))*dx
    

# define the initial solution
for i in range(nx_cells):

    for j in range(ny_cells):
    
        dist = np.sqrt( (X_cent[j,i]-7.0)**2 + (Y_cent[j,i]-0.0)**2 )

        Z_cent[j,i] = 0.25 * ( Z[j,i] + Z[j+1,i] + Z[j,i+1] + Z[j+1,i+1] )  

        if ( dist <= 1.85 ):
    
            W_cent[j,i] = np.sqrt( 1.85**2 - dist**2 ) + Z_cent[j,i]
            U_cent[j,i] = 0.0
            V_cent[j,i] = 0.0
    
        else: 
 
            W_cent[j,i] = Z_cent[j,i]
            U_cent[j,i] = 0.0
            V_cent[j,i] = 0.0

        H_cent[j,i] = W_cent[j,i] - Z_cent[j,i]


# create a figure for the plot
fig = plt.figure()
ax = fig.gca(projection='3d')
#plt.ylim([-0.1,1.5])
#plt.xlim([-0.25,1.75])

# plot the initial solution and call "line" the plot
# surf = ax.plot_surface(X,Y,Z)
# surf2 = ax.plot_surface(X_cent,Y_cent,W_cent)

X_cent1d = X_cent.flatten()
Y_cent1d = Y_cent.flatten()
H_cent1d = H_cent.flatten()
W_cent1d = W_cent.flatten()

idx = np.argwhere(H_cent1d>0)

idx = np.union1d(idx,idx+1)
idx = np.union1d(idx,idx-1)
idx = np.union1d(idx,idx+nx_cells)
idx = np.union1d(idx,idx-nx_cells)



# Triangulate parameter space to determine the triangles
tri = mtri.Triangulation(X_cent1d[idx].flatten(), Y_cent1d[idx].flatten())

# Plot the surface.  The triangles in parameter space determine which x, y, z
# points are connected by an edge.
ax.plot_trisurf(X_cent1d[idx].flatten(), Y_cent1d[idx].flatten(), W_cent1d[idx].flatten(), triangles=tri.triangles,edgecolor='none')

ax.plot_surface(X_cent, Y_cent, W_cent, alpha=0.3)
cset = ax.contour(X_cent, Y_cent, H_cent, zdir='z', offset=-5, cmap=cm.coolwarm)
cset = ax.contour(X_cent, Y_cent, H_cent, zdir='x', offset=0, cmap=cm.coolwarm)
cset = ax.contour(X_cent, Y_cent, H_cent, zdir='y', offset=7, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(0, 30)
ax.set_ylabel('Y')
ax.set_ylim(-15, 15)
ax.set_zlabel('Z')
ax.set_zlim(-5, 10)
#ax.axis('equal')

# fig.tight_layout()
plt.show()    


# create topography file
header = "ncols     %s\n" % nx_points
header += "nrows    %s\n" % ny_points
header += "xllcorner " + str(x_min-0.5*dx) +"\n"
header += "yllcorner " + str(y_min-0.5*dx) +"\n"
header += "cellsize " + str(dx) +"\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, Z, header=header, fmt='%1.12f',comments='')

# create intial solution file
init_file = 'example_2D_0000.q_2d'


for i in range(ny_cells):

    q0 = np.zeros((6+n_solid,nx_cells))

    q0[0,:] = X_cent[i,:]
    q0[1,:] = Y_cent[i,:]
    q0[2,:] = rho_m * H_cent[i,:]
    q0[3,:] = 0.0
    q0[4,:] = 0.0
    q0[5,:] = rho_m * H_cent[i,:] * SP_HEAT_MIX * T 

    for j in range(n_solid):
        q0[6+j,:] = rho_s*H_cent[i,:]*alfas/n_solid


    if ( i==0 ):
        
        with open(init_file, "w+") as file:
            np.savetxt(file, np.transpose(q0), fmt='%19.12e') 

    else:

        with open(init_file, "a") as file:
            np.savetxt(file, np.transpose(q0), fmt='%19.12e') 


    with open(init_file,'a') as file:
        file.write(' \n')



# Read in the file
with open('SW_VAR_DENS_MODEL.template', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'example2D')
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_min', str(x_min))
filedata = filedata.replace('y_min', str(y_min))
filedata = filedata.replace('nx_cells', str(nx_cells))
filedata = filedata.replace('ny_cells', str(ny_cells))
filedata = filedata.replace('dx', str(dx))


# Write the file out again
with open('SW_VAR_DENS_MODEL.inp', 'w') as file:
  file.write(filedata)

