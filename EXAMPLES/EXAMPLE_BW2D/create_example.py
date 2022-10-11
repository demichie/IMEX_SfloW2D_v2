#!/usr/bin/env python
"""
% This function 

"""
import numpy as np
import time
import sys

if len(sys.argv)==9: 

    a = sys.argv[1]

    if a.isdigit():

        nx_cells = int(a)
        ny_cells = nx_cells
 
    else:
 
        sys.exit()

    a = sys.argv[2]
    try:
        half_domain_size = float(a)
    except ValueError:
        print("You must enter a float for domain size "+a)

    a = sys.argv[3]
    try:
        r_init = float(a)
    except ValueError:
        print("You must enter a float for source radius "+a)

    a = sys.argv[4]
    try:
        slope = float(a)
    except ValueError:
        print("You must enter a float for the slope "+a)

    a = sys.argv[5]
    try:
        T_init = float(a)
    except ValueError:
        print("You must enter a float for source temperature: "+a)

    a = sys.argv[6]
    try:
        vel_init = float(a)
    except ValueError:
        print("You must enter a float for source radial velocity: "+a)

    a = sys.argv[7]
    try:
        alfas_init = float(a)
    except ValueError:
        print("You must enter a float for solid fraction")

    a = sys.argv[8]
    try:
        particle_size = float(a)
    except ValueError:
        print("You must enter a float for particle_size")

else:

    print('Please provide five arguments:\n')
    print('1) Number of cells\n')
    print('2) Half-domain size (>0)\n')
    print('3) Source radius (>0)\n')
    print('4) Slope (degrees, >=0)\n')
    print('5) Temperature (>0)\n')
    print('6) Radial velocity (>0)\n')
    print('7) Solid volume fraction (0,1)\n')
    print('8) Particle size (>0)\n')
    sys.exit()


n_solid = 1

# Define the boundaries x_left and x_right of the spatial domain
x_min = -half_domain_size
x_max = half_domain_size

y_min = -half_domain_size
y_max = half_domain_size


# Define the number n_points of points of the grid
nx_points  = nx_cells+1

# Define the grid stepsize dx
dx = ( x_max - x_min ) / ( nx_cells )

# print('dx',dx,x_max - x_min, nx_cells)

# Define the array x of the grid points
x = np.linspace(x_min,x_max,nx_points)

x_cent = np.linspace(x_min+0.5*dx,x_max-0.5*dx,nx_cells)
x_stag = np.linspace(x_min,x_max,nx_points)

dy = dx

ny_points = ny_cells+1

n_cells = nx_cells * ny_cells

# Define the array x of the grid points
y = np.linspace(y_min,y_max,ny_points)

y_cent = np.linspace(y_min+0.5*dy,y_max-0.5*dy,ny_cells)
y_stag = np.linspace(y_min,y_max,ny_points)

X, Y = np.meshgrid(x, y)
X_cent, Y_cent = np.meshgrid(x_cent, y_cent)
X_vrtx , Y_vrtx = np.meshgrid(x_stag,y_stag)
# print X.shape
# print X_cent.shape

Z_cent = -np.arctan(slope/180.0*np.pi)*np.sqrt( (X_cent)**2 + (Y_cent)**2 )

Z_vrtx = -np.arctan(slope/180.0*np.pi)*np.sqrt( (X_vrtx)**2 + (Y_vrtx)**2 )

Z_vrtx -= np.amin(Z_vrtx) 

H_cent = np.zeros_like(X_cent)
U_cent = np.zeros_like(X_cent)
V_cent = np.zeros_like(X_cent)



# create topography file
header = "ncols     %s\n" % nx_points
header += "nrows    %s\n" % ny_points
header += "xllcorner " + str(x_min-0.5*dx) +"\n"
header += "yllcorner " + str(y_min-0.5*dx) +"\n"
header += "cellsize " + str(dx) +"\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, Z_vrtx, header=header, fmt='%1.12f',comments='')

# create intial solution file
init_file = 'example_bw2d_0000.q_2d'


for i in range(ny_cells):

    q0 = np.zeros((6+n_solid,nx_cells))

    q0[0,:] = X_cent[i,:]
    q0[1,:] = Y_cent[i,:]
    q0[2,:] = 0.0
    q0[3,:] = 0.0
    q0[4,:] = 0.0
    q0[5,:] = 0.0 

    for j in range(n_solid):
        q0[6+j,:] = 0.0


    if ( i==0 ):
        
        with open(init_file, "w+") as file:
            np.savetxt(file, np.transpose(q0), fmt='%19.12e') 

    else:

        with open(init_file, "a") as file:
            np.savetxt(file, np.transpose(q0), fmt='%19.12e') 


    with open(init_file,'a') as file:
        file.write(' \n')



# Read in the file
with open('IMEX_SfloW2D.template', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'exampleRS')
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_min', str(x_min)+'D0')
filedata = filedata.replace('y_min', str(y_min)+'D0')
filedata = filedata.replace('nx_cells', str(nx_cells))
filedata = filedata.replace('ny_cells', str(ny_cells))
filedata = filedata.replace('dx', str(dx)+'D0')
filedata = filedata.replace('source_r', str(r_init)+'D0')
filedata = filedata.replace('source_T', str(T_init)+'D0')
filedata = filedata.replace('source_vel', str(vel_init)+'D0')
filedata = filedata.replace('source_alphas', str(alfas_init))
filedata = filedata.replace('particle_size', str(particle_size))


# Write the file out again
with open('IMEX_SfloW2D.inp', 'w') as file:
  file.write(filedata)







