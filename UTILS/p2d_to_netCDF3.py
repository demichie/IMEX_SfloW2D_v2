#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import time
import sys
import os.path
import netCDF4 


if len(sys.argv)==2: 
 
    bakfile = sys.argv[1]

    os.path.isfile(bakfile) 

else:

    print('Please provide bak file name:\n')
    print('1) File name\n')
    sys.exit()


print(bakfile)

with open(bakfile) as fp:  

   for cnt, line in enumerate(fp):

       if "T_START" in line:
           t_start_str= line.replace('T_START=','')
           t_start_str= t_start_str.replace(',','')
           t_start = np.float(t_start_str)
           print("t_start",t_start)

       if "T_END" in line:
           t_end_str= line.replace('T_END=','')
           t_end_str= t_end_str.replace(',','')
           t_end = np.float(t_end_str)
           print("t_start",t_start)

       if "DT_OUTPUT" in line:
           dt_output_str = line.replace('DT_OUTPUT=','')
           dt_output_str = dt_output_str.replace(',','')
           dt_output = np.float(dt_output_str)
           print("dt_output",dt_output)

       if "COMP_CELLS_X" in line:
           comp_cells_x_str= line.replace('COMP_CELLS_X=','')
           comp_cells_x_str= comp_cells_x_str.replace(',','')
           nx = np.int(comp_cells_x_str)
           print("comp_cells_x",nx)

       if "COMP_CELLS_Y" in line:
           comp_cells_y_str= line.replace('COMP_CELLS_Y=','')
           comp_cells_y_str= comp_cells_y_str.replace(',','')
           ny = np.int(comp_cells_y_str)
           print("comp_cells_y",ny)

       if ("X0" in line) and not("RUNOUT" in line):
           x0_str= line.replace('X0=','')
           x0_str= x0_str.replace(',','')
           x0 = np.float(x0_str)
           print("x0",x0)

       if ("Y0" in line) and not("RUNOUT" in line):
           y0_str= line.replace('Y0=','')
           y0_str= y0_str.replace(',','')
           y0 = np.float(y0_str)
           print("y0",y0)

       if "N_SOLID" in line:
           nsolid_str= line.replace('N_SOLID=','')
           nsolid_str= nsolid_str.replace(',','')
           nsolid = np.int(nsolid_str)
           print("nsolid",nsolid)


n_output = np.int((t_end-t_start)/dt_output)

print('n_output',n_output+1)

filename_fix = bakfile.replace('.bak','')


ncfilename = filename_fix+'.nc'

ncfile = netCDF4.Dataset(ncfilename,mode='w',format='NETCDF4_CLASSIC') 

x_dim = ncfile.createDimension('x', nx) 
y_dim = ncfile.createDimension('y', ny) 
time_dim = ncfile.createDimension('time', None) # unlimited axis (can be appended to).


ncfile.title=filename_fix+' output'

ncfile.Conventions = "CF-1.0"
ncfile.subtitle="My model data subtitle"
ncfile.anything="write anything"


x = ncfile.createVariable('x', np.float64, ('x',))
x.long_name = 'x dim'
x.units = 'meters'

y = ncfile.createVariable('y', np.float64, ('y',))
y.long_name = 'y dim'
y.units = 'meters'

time = ncfile.createVariable('time', np.float64, ('time',))
time.long_name = 'Time'
time.units = 'seconds'

h = ncfile.createVariable('h',np.float64,('time','y','x')) # note: unlimited dimension is leftmost
h.standard_name = 'flow thickness' # this is a CF standard name
h.units = 'meters' # degrees Kelvin

b = ncfile.createVariable('b',np.float64,('time','y','x'),zlib=True) 
b.standard_name = 'topography' # this is a CF standard name
b.units = 'meters' # degrees Kelvin

ux = ncfile.createVariable('ux',np.float64,('time','y','x'),zlib=True) 
ux.standard_name = 'x-velocity' # 
ux.units = 'meters/second' # 

uy = ncfile.createVariable('uy',np.float64,('time','y','x'),zlib=True) 
uy.standard_name = 'y-velocity' # 
uy.units = 'meters/second' # 

T = ncfile.createVariable('T',np.float64,('time','y','x'),zlib=True) 
T.standard_name = 'temperature' # this is a CF standard name
T.units = 'kelvin' # degrees Kelvin

W = ncfile.createVariable('W',np.float64,('time','y','x'),zlib=True) 
W.standard_name = 'free surface' # 
W.units = 'meters' # 

for i in range(nsolid):
    globals()['dep'+'_{0:04}'.format(i)] = ncfile.createVariable('dep'+'_{0:04}'.format(i),np.float64,('time','y','x'),zlib=True) 
    globals()['dep'+'_{0:04}'.format(i)].standard_name = 'deposit' # 
    globals()['dep'+'_{0:04}'.format(i)].units = 'meters' # 

    globals()['alphas'+'_{0:04}'.format(i)] = ncfile.createVariable('alphas'+'_{0:04}'.format(i),np.float64,('time','y','x'),zlib=True) 
    globals()['alphas'+'_{0:04}'.format(i)].standard_name = 'solid volume fraction' # 
    globals()['alphas'+'_{0:04}'.format(i)].units = '' # 

rhom = ncfile.createVariable('rhom',np.float64,('time','y','x'),zlib=True) 
rhom.standard_name = 'flow density' # 
rhom.units = 'kilograms/second' # 


X = np.zeros((ny,nx))
Y = np.zeros((ny,nx))



for i_output in range(0,n_output+1):

    filename = filename_fix+'_{0:04}'.format(i_output)+'.p_2d'
    print(filename)

    data = np.loadtxt(filename,skiprows=0)

    time[i_output] = dt_output*(i_output)

    h[i_output,:,:] = data[:,2].reshape((ny,nx))
    ux[i_output,:,:] = data[:,3].reshape((ny,nx))
    uy[i_output,:,:] = data[:,4].reshape((ny,nx))
    b[i_output,:,:] = data[:,5].reshape((ny,nx))
    W[i_output,:,:] = data[:,6].reshape((ny,nx))
    for i in range(nsolid):
        globals()['alphas'+'_{0:04}'.format(i)][i_output,:,:] = data[:,7+i].reshape((ny,nx))
    T[i_output,:,:] = data[:,7+nsolid].reshape((ny,nx))
    rhom[i_output,:,:] = data[:,8+nsolid].reshape((ny,nx))
    for i in range(nsolid):
        globals()['dep'+'_{0:04}'.format(i)][i_output,:,:] = data[:,10+nsolid+i].reshape((ny,nx))


X[:,:] = data[:,0].reshape((ny,nx))
Y[:,:] = data[:,1].reshape((ny,nx))

x[:] = X[0,:]
y[:] = Y[:,0]

print(ncfile)

ncfile.close(); print('Dataset is closed!')



