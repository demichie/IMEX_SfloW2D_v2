#!/usr/bin/env python
"""
% This function

"""
import numpy as np
import sys

if len(sys.argv) == 3:

    print('Number of cells')
    a = sys.argv[1]

    if a.isdigit():

        n_cells = int(a)
        print(n_cells)

    else:

        sys.exit()

    a = sys.argv[2]
    if a == 'true':
        plot_flag = True
    elif a == 'false':
        plot_flag = False
    else:
        raise ValueError("You must enter true or false: " + a)

else:

    print('Please provide two arguments:\n')
    print('1) Number of cells\n')
    print('2) Plot flag (true or false)\n')
    sys.exit()

rho_l = 1000.0
SP_HEAT_L = 4200.0
T = 300.0

# Define the boundaries x_left and x_right of the spatial domain
x_left = 0.0
x_right = 25.0

# Define the number n_points of points of the grid
n_points = n_cells + 1

# Define the grid stepsize dx
dx = (x_right - x_left) / (n_points - 1)

# Define the array x of the grid points
x = np.linspace(x_left, x_right, n_points)

x_cent = np.linspace(x_left + 0.5 * dx, x_right - 0.5 * dx, n_cells)

B = np.zeros((n_points, 1))
w = np.zeros((n_points, 1))
u = np.zeros((n_points, 1))

B_cent = np.zeros_like(x_cent)
w_cent = np.zeros_like(x_cent)
u_cent = np.zeros_like(x_cent)

slope = 10.0
B_disc = 10.0

h_left = 0.66
h_right = 0.66

u_left = 0.0
u_right = 0.0

w_left = 2.0
w_right = 2.0

init_sol_disc = -10.0

# define the topography
for i in range(n_points):

    if (np.abs(x[i] - B_disc) <= 2.0):

        B[i] = 0.2 - 0.05 * (x[i] - 10.0)**2

    else:

        B[i] = 0.0

# define the topography
for i in range(n_points):

    if x[i] < init_sol_disc:

        w[i] = w_left
        u[i] = u_left

    elif x[i] == init_sol_disc:

        w[i] = 0.5 * (w_left + w_right)
        u[i] = 0.5 * (u_left + u_right)

    else:

        w[i] = w_right
        u[i] = u_right

# define the initial solution
for i in range(n_cells):

    B_cent[i] = 0.5 * (B[i] + B[i + 1])
    w_cent[i] = 0.5 * (w[i] + w[i + 1])
    u_cent[i] = 0.5 * (u[i] + u[i + 1])

# create topography file
header = "ncols     %s\n" % n_points
header += "nrows    %s\n" % 1
header += "xllcorner " + str(x_left - 0.5 * dx) + "\n"
header += "yllcorner " + str(0 - 0.5 * dx) + "\n"
header += "cellsize " + str(dx) + "\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full,
           np.transpose(B),
           header=header,
           fmt='%1.12f',
           comments='')

# create intial solution file
q0 = np.zeros((6, n_cells))

q0[0, :] = x_cent
q0[1, :] = 0.0
q0[2, :] = rho_l * (w_cent - B_cent)
q0[3, :] = rho_l * (w_cent - B_cent) * u_cent
q0[4, :] = 0.0
q0[5, :] = rho_l * (w_cent - B_cent) * (SP_HEAT_L * T + 0.5 * u_cent**2)

init_file = 'exampleSubcritical_0000.q_2d'

np.savetxt(init_file, np.transpose(q0), fmt='%19.12e')

with open(init_file, 'a') as file:
    file.write('\n')

# Read in the file
with open('IMEX_SfloW2D.template', 'r') as file:
    filedata = file.read()

runname = 'exampleSub_' + str(n_cells)

# Replace the target string
filedata = filedata.replace('runname', runname)
filedata = filedata.replace('restartfile', init_file)
filedata = filedata.replace('x_left', str(x_left))
filedata = filedata.replace('n_cells', str(n_cells))
filedata = filedata.replace('dx', str(dx))

# Write the file out again
with open('IMEX_SfloW2D.inp', 'w') as file:
    file.write(filedata)

if (plot_flag):

    import matplotlib.pyplot as plt

    # create a figure for the plot
    fig, ax = plt.subplots()
    plt.ylim([-0.1, 5.1])
    plt.xlim([x_left, x_right])

    # plot the initial solution and call "line" the plot
    line1, = ax.plot(x, B)
    line2, = ax.plot(x_cent, w_cent, '-g')

    plt.show()
