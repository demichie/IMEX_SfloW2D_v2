#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
from mpl_toolkits.mplot3d import Axes3D                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys
import os.path
import matplotlib.tri as mtri
from matplotlib import cm


if len(sys.argv)==2: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 

else:

    print('Please provide two or three arguments:\n')
    print('1) File name\n')
    print('2) 1st Variable to plot: h,hB,B,u,v\n')
    print('2) 2nd Variable to plot: h,hB,B,u,v\n')
    sys.exit()

data = np.loadtxt(filename,skiprows=0)

x = data[:,0]
y = data[:,1]
h = data[:,2]
u = data[:,3]
v = data[:,4]
b = data[:,5]
w = data[:,6]

mag_vel = np.sqrt(u**2+v**2)

x0_idx = np.asarray((np.where(data[:,0]==data[0,0])))

ny_cells = x0_idx[0,1]
nx_cells = data.shape[0] / ny_cells

nx_cells = nx_cells.astype(int)
ny_cells = ny_cells.astype(int)

X_cent = x.reshape((nx_cells,ny_cells))
Y_cent = y.reshape((nx_cells,ny_cells))
H_cent = h.reshape((nx_cells,ny_cells))
B_cent = b.reshape((nx_cells,ny_cells))
W_cent = w.reshape((nx_cells,ny_cells))
MAG_VEL_cent = mag_vel.reshape((nx_cells,ny_cells))

# create a figure for the plot
fig = plt.figure()
ax = fig.gca(projection='3d')

X_cent1d = X_cent.flatten()
Y_cent1d = Y_cent.flatten()
H_cent1d = H_cent.flatten()
W_cent1d = W_cent.flatten()

idx = np.argwhere(H_cent1d>1e-3)

# Triangulate parameter space to determine the triangles
tri = mtri.Triangulation(X_cent1d[idx].flatten(), Y_cent1d[idx].flatten())

# Plot the surface.  The triangles in parameter space determine which x, y, z
# points are connected by an edge.
ax.plot_trisurf(X_cent1d[idx].flatten(), Y_cent1d[idx].flatten(), W_cent1d[idx].flatten(), triangles=tri.triangles,edgecolor='none')

ax.plot_surface(X_cent, Y_cent, W_cent, alpha=0.2)

#cset = ax.contour(X_cent, Y_cent, H_cent, zdir='z', offset=-5, cmap=cm.coolwarm)
cset = ax.contour(X_cent, Y_cent, MAG_VEL_cent, zdir='z', offset=-5, cmap=cm.coolwarm)

#cset = ax.contour(X_cent, Y_cent, H_cent, zdir='x', offset=0, cmap=cm.coolwarm)
#cset = ax.contour(X_cent, Y_cent, H_cent, zdir='y', offset=7, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_xlim(0, 30)
ax.set_ylabel('Y')
ax.set_ylim(-15, 15)
ax.set_zlabel('Z')
ax.set_zlim(-5, 10)
#ax.axis('equal')

labels = ['line1', 'line2','line3','line4',
           'line5', 'line6']

levels = cset.levels

for i in range(len(levels)):
    cset.collections[i].set_label(str(levels[i]))

plt.legend(loc='upper left',title='mag(u)')
# fig.tight_layout()
plt.show()    


