#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys

if len(sys.argv)==8: 

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        n_cells = int(a)
        print n_cells
 
    else:
 
        sys.exit()

    a = sys.argv[4]
    try:
        alfas = float(a)
    except ValueError:
        print("You must enter a float")

    a = sys.argv[5]
    try:
        erosion_coeff = float(a)
    except ValueError:
        print("You must enter a float for erosion coefficient: "+a)

    a = sys.argv[6]
    try:
        settling_vel = float(a)
    except ValueError:
        print("You must enter a float for settling velocity: "+a)


    a = sys.argv[7]
    try:
        T = float(a)
    except ValueError:
        print("You must enter a float for temperature: "+a)


else:

    print('Please provide three arguments:\n')
    print('1) Number of cells\n')
    print('2) Variables to reconstruct: phys or cons\n')
    print('3) Order of the RK scheme\n')
    print('4) Solid volume fraction (0,1)\n')
    print('5) Erosion coefficient (>0)\n')
    print('6) Settling velocity (>0)\n')
    print('7) Temperature (>0)\n')
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
x_left = 0.0
x_right = 50000.0


# Define the number n_points of points of the grid
n_points  = n_cells+1

# Define the grid stepsize dx
dx = ( x_right - x_left ) / ( n_points - 1 )


# Define the array x of the grid points
x = np.linspace(x_left,x_right,n_points)

x_cent = np.linspace(x_left+0.5*dx,x_right-0.5*dx,n_cells)

B = np.zeros((n_points,1))
w = np.zeros((n_points,1))
u = np.zeros((n_points,1))

B_cent = np.zeros_like(x_cent)
w_cent = np.zeros_like(x_cent)
u_cent = np.zeros_like(x_cent)

slope = 0.0
B_disc = 0.0


h_left = 200.0
h_right = 0.0

u_left = 0.0
u_right = 0.0

init_sol_disc = 1000.0

# define the topography
for i in range(n_points):

    B[i] = np.abs(x[i]) * np.sin(np.deg2rad(slope))
    

B = B - np.amin(B)

# define the topography
for i in range(n_points):

    if x[i]<init_sol_disc:
    
        w[i] = B[i]+h_left
        u[i] = u_left
    
    elif x[i]==init_sol_disc:
    
        w[i] = B[i]+0.5*(h_left+h_right)
        u[i] = 0.5*(u_left+u_right)

    else:
    
        w[i] = B[i]+h_right
        u[i] = u_right

 
# define the initial solution
for i in range(n_cells):

    B_cent[i] = 0.5*(B[i]+B[i+1])
    w_cent[i] = 0.5*(w[i]+w[i+1])
    u_cent[i] = 0.5*(u[i]+u[i+1])


# create topography file
header = "ncols     %s\n" % n_points
header += "nrows    %s\n" % 1
header += "xllcorner " + str(x_left-0.5*dx) +"\n"
header += "yllcorner " + str(0-0.5*dx) +"\n"
header += "cellsize " + str(dx) +"\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, np.transpose(B), header=header, fmt='%1.12f',comments='')

# create intial solution file
q0 = np.zeros((6+n_solid,n_cells))

q0[0,:] = x_cent
q0[1,:] = 0.0
q0[2,:] = rho_m * (w_cent-B_cent)
q0[3,:] = rho_m * (w_cent-B_cent)*u_cent
q0[4,:] = 0.0
q0[5,:] = rho_m * (w_cent-B_cent)* ( SP_HEAT_MIX * T + 0.5 * u_cent**2 )

for i in range(n_solid):
    q0[6+i,:] = rho_s*(w_cent-B_cent)*alfas/n_solid

init_file = 'example_BOX_0000.q_2d'

np.savetxt(init_file, np.transpose(q0), fmt='%19.12e') 

with open(init_file,'a') as file:
    file.write('\n')

# Read in the file
with open('IMEX_SfloW2D.template', 'r') as file :
  filedata = file.read()

runname= 'exampleBOX_'+str(n_cells)

# Replace the target string
filedata = filedata.replace('runname', runname)
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_left', str(x_left))
filedata = filedata.replace('n_cells', str(n_cells))
filedata = filedata.replace('dx', str(dx))
filedata = filedata.replace('erosionCoeff', str(erosion_coeff))
filedata = filedata.replace('settlingVel', str(settling_vel))

if sys.argv[2]=='cons':

    print('Linear reconstruction of conservative variables (h+B,hu,hv)')
    filedata = filedata.replace('bc2', 'HU')
    filedata = filedata.replace('recvar', 'cons')

else:

    print('Linear reconstruction of physical variables (h+B,u,v)')
    filedata = filedata.replace('bc2', 'U')
    filedata = filedata.replace('recvar', 'phys')

filedata = filedata.replace('order', sys.argv[3])


# Write the file out again
with open('IMEX_SfloW2D.inp', 'w') as file:
  file.write(filedata)

# create a figure for the plot
fig, ax = plt.subplots()
#plt.ylim([-0.1,5.1])
plt.xlim([x_left,x_right])

# plot the initial solution and call "line" the plot
line1, = ax.plot(x,B)
line2, = ax.plot(x_cent,w_cent,'-g')


plt.show()
