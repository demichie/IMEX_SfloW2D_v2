#!/usr/bin/env python
"""
% This function 

"""
import numpy as np
import time
import sys

if len(sys.argv)==7: 

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        nx_cells = int(a)
 
    else:
 
        sys.exit()

    a = sys.argv[2]
    try:
        r_init = float(a)
    except ValueError:
        print("You must enter a float for source radius "+a)

    a = sys.argv[3]
    try:
        h_init = float(a)
    except ValueError:
        print("You must enter a float for initial thickness: "+a)

    a = sys.argv[4]
    try:
        T_init = float(a)
    except ValueError:
        print("You must enter a float for source temperature: "+a)

    a = sys.argv[5]
    try:
        Ri_init = float(a)
    except ValueError:
        print("You must enter a float for Richardson number: "+a)

    a = sys.argv[6]
    try:
        xs_init = float(a)
    except ValueError:
        print("You must enter a float for solid fraction")

else:

    print('Please provide three arguments:\n')
    print('1) Number of cells\n')
    print('2) Source radius (>0)\n')
    print('3) Initial thickness (>0)\n')
    print('4) Temperature (>0)\n')
    print('5) Radial velocity (>0)\n')
    print('6) Solid mass fraction (0,1)\n')
    sys.exit()


n_solid = 1
rho_s = 2000.0
pres = 101300.0
sp_gas_const = 287.051
rho_g =  pres / ( sp_gas_const * T_init )
xg_init = 1.0 - xs_init

alfag_init = xg_init * rho_s / ( rho_g - xg_init * ( rho_g - rho_s ) )
alfas_init = 1.0 - alfag_init

rho_mix = 1.0/(xs_init/rho_s + (1.0-xs_init)/rho_g)

print('alfas_init',alfas_init)
print('rho_mix',alfas_init*rho_s + (1.0-alfas_init)*rho_g)
print('rho_mix',rho_mix)

T_ambient = 300.0
rho_atm = pres / ( sp_gas_const * T_ambient )
print('rho_atm',rho_atm)
grav = 9.81

vel_init = np.sqrt( grav * h_init * ( rho_mix - rho_atm ) / ( rho_mix * Ri_init ) )

print('vel_init',vel_init)

# Define the boundaries x_left and x_right of the spatial domain
x_min = -10000.0
x_max = 10000.0

y_min = -10000.0
y_max = 10000.0


# Define the number n_points of points of the grid
nx_points  = nx_cells+1

# Define the grid stepsize dx
dx = ( x_max - x_min ) / ( nx_cells )

# print('dx',dx,x_max - x_min, nx_cells)

# Define the array x of the grid points
x = np.linspace(x_min,x_max,nx_points)

x_cent = np.linspace(x_min+0.5*dx,x_max-0.5*dx,nx_cells)

x_half = 0.5 * ( x_min + x_max)
y_half = 0.5 * ( y_min + y_max)

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

    Z[:,i] = 1.0
    
# define the initial solution
for i in range(nx_cells):

    for j in range(ny_cells):
    
        dist = np.sqrt( (X_cent[j,i]-7.0)**2 + (Y_cent[j,i]-0.0)**2 )

        Z_cent[j,i] = 0.25 * ( Z[j,i] + Z[j+1,i] + Z[j,i+1] + Z[j+1,i+1] )  

        if ( dist <= 1.85 ):
    
            W_cent[j,i] = Z_cent[j,i]
            U_cent[j,i] = 0.0
            V_cent[j,i] = 0.0
    
        else: 
 
            W_cent[j,i] = Z_cent[j,i]
            U_cent[j,i] = 0.0
            V_cent[j,i] = 0.0

        H_cent[j,i] = W_cent[j,i] - Z_cent[j,i]




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
init_file = 'example_RS_0000.q_2d'


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
with open('SW_VAR_DENS_MODEL.template', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'exampleRS')
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_min', str(x_min)+'D0')
filedata = filedata.replace('y_min', str(y_min)+'D0')
filedata = filedata.replace('nx_cells', str(nx_cells))
filedata = filedata.replace('ny_cells', str(ny_cells))
filedata = filedata.replace('dx', str(dx)+'D0')
filedata = filedata.replace('source_x', str(x_half)+'D0')
filedata = filedata.replace('source_y', str(y_half)+'D0')
filedata = filedata.replace('source_h', str(h_init)+'D0')
filedata = filedata.replace('source_r', str(r_init)+'D0')
filedata = filedata.replace('source_T', str(T_init)+'D0')
filedata = filedata.replace('source_vel', str(vel_init)+'D0')
filedata = filedata.replace('source_alphas', str(alfas_init))


# Write the file out again
with open('SW_VAR_DENS_MODEL.inp', 'w') as file:
  file.write(filedata)







