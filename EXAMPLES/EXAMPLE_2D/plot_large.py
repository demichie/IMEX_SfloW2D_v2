#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import time
import sys
import os.path

print len(sys.argv)

if len(sys.argv)==3: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 

    var = sys.argv[2]

else:

    print('Please provide two or three arguments:\n')
    print('1) File name\n')
    print('2) 1st Variable to plot: u,T,rho_m\n')
    sys.exit()

data = np.loadtxt(filename,skiprows=0)



print data.shape

x = data[:,0]
y = data[:,1]
h = data[:,2]
u = data[:,3]
v = data[:,4]
b = data[:,5]
w = data[:,6]
alphas = data[:,7]
T = data[:,8]
rho_m = data[:,9]
red_grav = data[:,10]
mag_vel = np.sqrt(u**2+v**2)

if var=='u':

    z2d = mag_vel

elif var=='T':

    z2d = T

elif var=='alphas':

    z2d = alphas

elif var=='rho_m':

    z2d = rho_m

else:

    print('Please specify the variable to plot as 2nd argument: u,T,rho_m')
    sys.exit()


x0_idx = np.asarray((np.where(data[:,0]==data[0,0])))

ny_cells = x0_idx[0,1]
nx_cells = data.shape[0] / ny_cells

X_cent = x.reshape((nx_cells,ny_cells))
Y_cent = y.reshape((nx_cells,ny_cells))
H_cent = h.reshape((nx_cells,ny_cells))
B_cent = b.reshape((nx_cells,ny_cells))
W_cent = w.reshape((nx_cells,ny_cells))

W_cent = B_cent + H_cent

Wmin = np.min(W_cent)
Wmax = np.max(W_cent)

Wlin = np.around(np.linspace(Wmin,Wmax,num=5),decimals=2)

Wlin_str = []
for i in Wlin:
    Wlin_str.append(str(i))

W_cent = (W_cent-Wmin)/(Wmax-Wmin)

Wlin = np.around(np.linspace(0.0,1.0,num=5),decimals=2)


z2d_cent = z2d.reshape((nx_cells,ny_cells))

z2d_min = np.min(z2d_cent)
z2d_max = np.max(z2d_cent)

z2d_lin = np.around(np.linspace(z2d_min,z2d_max,num=5),decimals=4)
z2d_lin_str = []
for i in z2d_lin:
    z2d_lin_str.append(str(i))


z2d_cent = 0.001 * ( z2d_cent - z2d_min ) / ( z2d_max - z2d_min )

z2d_min = np.min(z2d_cent)
z2d_max = np.max(z2d_cent)


z2d_lin = np.around(np.linspace(z2d_min,z2d_max,num=5),decimals=7)


idx = np.ma.masked_where(H_cent<=0.0001,z2d_cent)
z2d_cent[np.where(np.ma.getmask(idx)==True)] = np.nan

import plotly
import plotly.graph_objects as go

import pandas as pd

trace2 = go.Surface(x=X_cent,y=Y_cent,z=z2d_cent, visible=True,opacity=1.0)

data = [trace2]

layout = go.Layout(title='Example2D', autosize=False,
                  width=1200, height=800,
                  margin=dict(l=65, r=50, b=65, t=90))


fig = go.Figure(data=data, 
        layout=layout)

fig.update_traces(colorbar=dict(
        title=var,
        titleside="top",
        tickmode="array",
        tickvals=z2d_lin,
        ticktext=z2d_lin_str,
        ticks="outside"
        ))

fig.add_surface(x=X_cent,y=Y_cent,z=W_cent,surfacecolor=H_cent, opacity=0.95,showscale=False,
              colorscale=[[0.0, "rgb(214,214,0)"],
                [0.005, "rgb(215,148,39)"],
                [1.0, "rgb(215,148,39)"]])

zratio = (Wmax-Wmin) / ( np.max(X_cent)-np.min(X_cent) )
yratio = ( np.max(Y_cent)-np.min(Y_cent) ) / ( np.max(X_cent)-np.min(X_cent) )

fig.update_layout(scene = dict(
                  zaxis = dict(
                  tickmode = 'array',
                  tickvals = Wlin,
                  ticktext = Wlin_str),
                  aspectmode = "manual",
                  aspectratio =dict(x=1,y=yratio,z=zratio))
                  )

camera = dict(
    eye=dict(x=1.0, y=-1.0, z=0.5)
)

fig.update_layout(scene_camera=camera)

fig.show()



