#!/usr/bin/env python
"""
% This function

"""
import numpy as np
import sys

if len(sys.argv) == 7:

    a = sys.argv[1]

    if a.isdigit():

        nx_cells = int(a)

    else:

        sys.exit()

    a = sys.argv[2]
    try:
        dy_init = float(a)
    except ValueError:
        print("You must enter a float for source extent: " + a)

    a = sys.argv[3]
    try:
        mfr_init = float(a)
    except ValueError:
        print("You must enter a float for source mass flow rate " + a)

    a = sys.argv[4]
    try:
        T_init = float(a)
    except ValueError:
        print("You must enter a float for source temperature: " + a)

    a = sys.argv[5]
    try:
        Ri_init = float(a)
    except ValueError:
        print("You must enter a float for Richardson number: " + a)

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
y1_init = -0.5 * dy_init
y2_init = 0.5 * dy_init

# Define the boundaries x_left and x_right of the spatial domain
x_min = -10000.0
x_max = 10000.0

y_min = -10000.0
y_max = 10000.0

x_min = -4000.0
x_max = 4000.0

y_min = -4000.0
y_max = 4000.0

# Define the number n_points of points of the grid
nx_points = nx_cells + 1

# Define the grid stepsize dx
dx = (x_max - x_min) / (nx_cells)

# print('dx',dx,x_max - x_min, nx_cells)

# Define the array x of the grid points
x = np.linspace(x_min, x_max, nx_points)

x_cent = np.linspace(x_min + 0.5 * dx, x_max - 0.5 * dx, nx_cells)

x_half = 0.5 * (x_min + x_max)
y_half = 0.5 * (y_min + y_max)

dy = dx

# print('dy',dy)
ny_half_cells = int(np.ceil(y_max / dy))
ny_cells = 2 * ny_half_cells
ny_points = ny_cells + 1

y_min = -dy * ny_half_cells
y_max = -y_min

print('Number of cells in the y-direction')
print(ny_cells)

n_cells = nx_cells * ny_cells

# print(y_min)
# print(y_max)

# Define the array x of the grid points
y = np.linspace(y_min, y_max, ny_points)

y_cent = np.linspace(y_min + 0.5 * dy, y_max - 0.5 * dy, ny_cells)

X, Y = np.meshgrid(x, y)
X_cent, Y_cent = np.meshgrid(x_cent, y_cent)

Z = np.zeros_like(X)

# define the topography
for i in range(nx_points - 1, -1, -1):

    Z[:, i] = 1.0

# Create topography file
header = "ncols     %s\n" % nx_points
header += "nrows    %s\n" % ny_points
header += "xllcorner " + str(x_min - 0.5 * dx) + "\n"
header += "yllcorner " + str(y_min - 0.5 * dx) + "\n"
header += "cellsize " + str(dx) + "\n"
header += "NODATA_value -9999\n"

output_full = 'topography_dem.asc'

np.savetxt(output_full, Z, header=header, fmt='%1.12f', comments='')

# Read in the template file
with open('IMEX_SfloW2D.template', 'r') as file:
    filedata = file.read()

# Replace the target string
filedata = filedata.replace('runname', 'example_SE2016')
filedata = filedata.replace('x_min', str(x_min) + 'D0')
filedata = filedata.replace('y_min', str(y_min) + 'D0')
filedata = filedata.replace('nx_cells', str(nx_cells))
filedata = filedata.replace('ny_cells', str(ny_cells))
filedata = filedata.replace('dx', str(dx) + 'D0')
filedata = filedata.replace('source_y1', str(y1_init) + 'D0')
filedata = filedata.replace('source_y2', str(y2_init) + 'D0')
filedata = filedata.replace('source_mfr', str(mfr_init) + 'D0')
filedata = filedata.replace('source_Ri', str(Ri_init) + 'D0')
filedata = filedata.replace('source_T', str(T_init) + 'D0')
filedata = filedata.replace('source_xs', str(xs_init) + 'D0')

# Write the input file
with open('IMEX_SfloW2D.inp', 'w') as file:
    file.write(filedata)
