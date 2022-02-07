import matplotlib.pyplot as plt
import numpy as np
from linecache import getline
import matplotlib.patches as mpatches
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker
import csv
import glob, os
from os.path import exists
from matplotlib.colors import LightSource

source1 = '../DEM/vulcano_1m_hillshade.asc'

        

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

for file in glob.glob("*.bak"):
    # print(file)
    bakfile = file

with open(bakfile) as fp:  

   for cnt, line in enumerate(fp):
   
       if "TOPOGRAPHY_FILE" in line:
       
           topography_file = line.replace('TOPOGRAPHY_FILE','')
           topography_file = topography_file.replace('=','')
           topography_file = topography_file.replace('"','')
           topography_file = topography_file.replace(',','')
           topography_file = topography_file.replace(' ','')
           topography_file = os.getcwd()+'/'+topography_file.rstrip()
           
       if "T_START" in line:
           t_start_str= line.replace('T_START','')
           t_start_str= t_start_str.replace('=','')
           t_start_str= t_start_str.replace(',','')
           t_start = float(t_start_str)
           # print("t_start =",t_start)

       if "T_END" in line:
           t_end_str= line.replace('T_END','')
           t_end_str= t_end_str.replace('=','')
           t_end_str= t_end_str.replace(',','')
           t_end = float(t_end_str)
           # print("t_end =",t_end)

       if "DT_OUTPUT" in line:
           dt_output_str = line.replace('DT_OUTPUT','')
           dt_output_str = dt_output_str.replace('=','')
           dt_output_str = dt_output_str.replace(',','')
           dt_output = float(dt_output_str)
           # print("dt_output =",dt_output)
           
       if "X0" in line:
           x0_str = line.replace('X0','')
           x0_str = x0_str.replace('=','')
           x0_str = x0_str.replace(',','')
           x0 = float(x0_str)
           print("x0 =",x0)    

       if "Y0" in line:
           y0_str = line.replace('Y0','')
           y0_str = y0_str.replace('=','')
           y0_str = y0_str.replace(',','')
           y0 = float(y0_str)
           print("y0 =",y0)    
           
       if "COMP_CELLS_X" in line:
           comp_cells_x_str = line.replace('COMP_CELLS_X','')
           comp_cells_x_str = comp_cells_x_str.replace('=','')
           comp_cells_x_str = comp_cells_x_str.replace(',','')
           comp_cells_x = float(comp_cells_x_str)
           print("comp_cells_x =",comp_cells_x)    

       if "COMP_CELLS_Y" in line:
           comp_cells_y_str = line.replace('COMP_CELLS_Y','')
           comp_cells_y_str = comp_cells_y_str.replace('=','')
           comp_cells_y_str = comp_cells_y_str.replace(',','')
           comp_cells_y = float(comp_cells_y_str)
           print("comp_cells_y =",comp_cells_y)    

       if "CELL_SIZE" in line:
           cell_size_str = line.replace('CELL_SIZE','')
           cell_size_str = cell_size_str.replace('=','')
           cell_size_str = cell_size_str.replace(',','')
           cell_size = float(cell_size_str)
           print("cell_size =",cell_size)    


plt_xmin = x0
plt_xmax = x0+comp_cells_x*cell_size
plt_ymin = y0
plt_ymax = y0+comp_cells_y*cell_size


# Parse the header using a loop and
# the built-in linecache module
hdr = [getline(topography_file, i) for i in range(1,7)]
values = [float(h.split(" ")[-1].strip()) \
 for h in hdr]
cols,rows,lx,ly,cell,nd = values
xres = cell
yres = cell * -1

# Load the dem into a numpy array
h = np.loadtxt(topography_file, skiprows=7)

# make these smaller to increase the resolution
dx, dy = 0.05, 0.05

x = np.arange(-3.0, 3.0, dx)
y = np.arange(-3.0, 3.0, dy)
X, Y = np.meshgrid(x, y)

# when layering multiple images, the images need to have the same
# extent.  This does not mean they need to have the same shape, but
# they both need to render to the same coordinate system determined by
# xmin, xmax, ymin, ymax.  Note if you use different interpolations
# for the images their apparent extent could be different due to
# interpolation edge effects



fig = plt.figure()
ax = fig.add_subplot(111)

ls = LightSource(azdeg=315, altdeg=45)

extent1 = [lx, lx+cols*cell, ly, ly+rows*cell]

im1 = plt.imshow(ls.hillshade(np.flipud(h),vert_exag=1.0,
                     dx=cell,dy=cell),cmap='gray',extent=extent1)

plt.xlim([plt_xmin, plt_xmax])
plt.ylim([plt_ymin, plt_ymax])

ax.set_aspect('equal', 'box')

# plt.show()


VTfile = bakfile.replace('.bak','_VT.txt')

print('VTfile',VTfile,os.path.exists(VTfile))
# print(VTfile)
if os.path.exists(VTfile):

    idx = []
    thickness =[]
    dyn_pressure = []

    with open(VTfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ', 
                        skipinitialspace=True)

        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                idx.append(row[0])
                thickness.append(row[1])
                dyn_pressure.append(row[2])
                line_count += 1
        # print(f'Processed {line_count} lines.')


    thickness_unique = np.unique(thickness)
    # print(thickness_unique)
    dyn_pressure_unique = np.sort(np.unique([float(i) for i in dyn_pressure]))
    # print(dyn_pressure_unique)





os.chdir(".")

counter = 0

min_arr_T = 1.E10
max_arr_T = -1.E10

min_arr_dep = 1.E10
max_arr_dep = -1.E10

min_arr_ers = 1.E10
max_arr_ers = -1.E10

min_arr_thk = 1.E10
max_arr_thk = -1.E10


for file in glob.glob("*.asc"):
    print(file)

    if "masked" in file:

        continue
        
    if "topography_dem" in file:

        continue
        
    if "area" in file:

        continue
        
    if "esri" in file:

        continue
        
    if "slope" in file:

        continue
        
    if "deposit" in file:

        continue
        
    if "_T_" in file:
    
        continue    
        
    elif "_dep_" in file:

        if "_0000." in file:
            continue

    elif "_ers_" in file:

        if "_0000." in file:
            continue

    if "_max" in file:

        t_output = t_end

    elif "_vt_" in file:

        t_output = t_end

    else:

        file_idx = int(file[file.index(".")-4:file.index(".")])
        t_output = t_start+file_idx*dt_output 

    t_string = 't='+str(t_output)+'s'
    print(t_string)


    if counter==0:

        source2 = file

        # Parse the header using a loop and
        # the built-in linecache module
        hdr = [getline(source2, i) for i in range(1,7)]
        values = [float(h.split(" ")[-1].strip()) \
         for h in hdr]
        cols,rows,lx,ly,cell,nd = values
        xres = cell
        yres = cell * -1

        # Load the dem into a numpy array
        arr = np.loadtxt(source2, skiprows=6)
        arr[arr==nd]=np.nan
        arr=arr.ravel()
        arr=np.flipud(arr.reshape((int(rows),int(cols))))
        
        if np.all(np.isnan(arr)):
        
            continue        

        if "_vt_" not in file:
            # compute the range of values to plot
            min_arr = np.nanmin(arr)
            max_arr = np.nanmax(arr)
        
            min_int = np.minimum(np.floor(np.log10(max_arr))-4,np.maximum(-5,np.floor(np.log10(min_arr))))
            max_int = np.minimum(1,np.floor(np.log10(max_arr)))

        if "_T_" in file:
        
            min_arr_T = np.minimum(min_arr_T,min_arr)
            max_arr_T = np.maximum(max_arr_T,max_arr)

            # print('temperature',min_arr_T,max_arr_T)

            min_int_T = np.minimum(np.floor(np.log10(max_arr_T))-4,np.maximum(-5,np.floor(np.log10(min_arr_T))))
            max_int_t = np.minimum(1,np.floor(np.log10(max_arr_T)))


        elif "_dep" in file:

            min_arr_dep = np.minimum(min_arr_dep,min_arr)
            max_arr_dep = np.maximum(max_arr_dep,max_arr)

            # print('deposit',min_arr_dep,max_arr_dep)

            min_int_dep = np.minimum(np.floor(np.log10(max_arr_dep))-4,np.maximum(-5,np.floor(np.log10(min_arr_dep))))
            max_int_dep = np.minimum(1,np.floor(np.log10(max_arr_dep)))

        elif "_ers" in file:

            min_arr_ers = np.minimum(min_arr_ers,min_arr)
            max_arr_ers = np.maximum(max_arr_ers,max_arr)

            # print('erosion',min_arr_ers,max_arr_ers)

            min_int_ers = np.minimum(np.floor(np.log10(max_arr_ers))-4,np.maximum(-5,np.floor(np.log10(min_arr_ers))))
            max_int_ers = np.minimum(1,np.floor(np.log10(max_arr_ers)))

        elif "_vt_" in file:
        
            print('')
        

        else:

            min_arr_thk = np.minimum(min_arr_thk,min_arr)
            max_arr_thk = np.maximum(max_arr_thk,max_arr)

            # print('thickness',min_arr_thk,max_arr_thk)

            min_int_thk = np.minimum(np.floor(np.log10(max_arr_thk))-4,np.maximum(-5,np.floor(np.log10(min_arr_thk))))
            max_int_thk = np.minimum(1,np.floor(np.log10(max_arr_thk)))
        

for file in glob.glob("*.asc"):
    # print('QUI',file)

    if "topography_dem" in file:

        continue
        
    if "masked" in file:

        continue
        
    if "area" in file:

        continue
        
    if "esri" in file:

        continue
        
    if "_T_" in file:

        continue
        
    if "slope" in file:

        continue
        
    if "deposit" in file:

        continue
        
    elif "_dep_" in file:

        if "_0000." in file:
            continue

    elif "_ers_" in file:

        if "_0000." in file:
            continue

    if "_max" in file:

        t_output = t_end

    elif "_vt_" in file:

        t_output = t_end

    else:

        file_idx = int(file[file.index(".")-4:file.index(".")])
        t_output = t_start+file_idx*dt_output
        
    t_string = 't={t_output:.1f}s'.format(t_output = t_output)
    # print('t_string',t_string,t_output)


    if counter==0:

        source2 = file
        print(source2)

        # Parse the header using a loop and
        # the built-in linecache module
        hdr = [getline(source2, i) for i in range(1,7)]
        values = [float(h.split(" ")[-1].strip()) \
         for h in hdr]
        cols,rows,lx,ly,cell,nd = values
        xres = cell
        yres = cell * -1

        # Load the dem into a numpy array
        arr = np.loadtxt(source2, skiprows=6)
        arr[arr==nd]=np.nan
        arr = np.flipud(arr)

        
        if np.all(np.isnan(arr)):
        
            continue        

        # when layering multiple images, the images need to have the same
        # extent.  This does not mean they need to have the same shape, but
        # they both need to render to the same coordinate system determined by
        # xmin, xmax, ymin, ymax.  Note if you use different interpolations
        # for the images their apparent extent could be different due to
        # interpolation edge effects

        x = np.linspace(lx,lx+(cols)*cell,int(cols+1))
        y = np.linspace(ly,ly+(rows)*cell,int(rows+1))
        
        X,Y = np.meshgrid(x,y)

        extent2 = lx, lx+cols*cell, ly, ly+rows*cell

        if "_vt_" not in file:
            # compute the range of values to plot
            min_arr = np.nanmin(arr)
            max_arr = np.nanmax(arr)
        
            min_int = np.minimum(np.floor(np.log10(max_arr))-4,np.maximum(-5,np.floor(np.log10(min_arr))))
            max_int = np.minimum(1,np.floor(np.log10(max_arr)))

        if "_T_" in file:

            # print('temperature')
            levels = np.linspace(min_arr_T,max_arr_T,10)
            label_str = 'Temperature (K)'

        elif "_dep" in file:

            # print('deposit')
            levels = 10.0**np.linspace(min_int,max_int_dep,2*int(max_int_dep-min_int_dep)+1)
            label_str = 'Deposit thickness (m)'

        elif "_ers" in file:

            # print('erosion')
            levels = 10.0**np.linspace(min_int,max_int_ers,2*int(max_int_ers-min_int_ers)+1)
            label_str = 'Erosion thickness (m)'

        elif "_vt_" in file:

            # print('dynamic pressure')
            levels = dyn_pressure_unique            
            lev = 0
            for dyn_lev in dyn_pressure_unique:

                arr[arr==lev]=dyn_lev
                lev +=1
            
            label_str = 'Dynamic pressure (Pa)'
            levels = np.append(levels,np.max(levels)*2)

        else:

            # print('thickness')
            levels = 10.0**np.linspace(min_int,max_int_thk,2*int(max_int_thk-min_int_thk)+1)
            label_str = 'Flow thickness (m)'

        
        cmap = plt.get_cmap('terrain_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        p1 = plt.pcolormesh(X,Y,arr, cmap=cmap, norm=norm,alpha=0.65)
        #p1 = plt.pcolormesh(x,y,arr, cmap=cmap, norm=norm,alpha=0.65)
        
        if "_vt_" in file:
            ticks  = 0.5*(levels[:-1]+levels[1:])
            clb = plt.colorbar(format=ticker.FuncFormatter(fmt),ticks=ticks)
            clb.set_ticklabels(np.char.mod('%d', levels[:-1]))
            
        else:
        
            clb = plt.colorbar(format=ticker.FuncFormatter(fmt))

        clb.set_label(label_str, labelpad=-40, y=1.05, rotation=0)

        txt_box = plt.text(0.1, 0.9,t_string,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes,bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 5})            
        
        plt.axis('scaled')
        plt.xlim([plt_xmin,plt_xmax])
        plt.ylim([plt_ymin,plt_ymax])
        plt.tight_layout()
        plt.savefig(source2.replace('asc','jpg'),dpi=200)

    else:

        source2 = file
        print(source2)
        
        # Load the dem into a numpy array
        arr = np.loadtxt(source2, skiprows=6)
        arr[arr==nd]=np.nan
        arr=np.flipud(arr)
        #arr=arr.ravel()
        #arr=np.flipud(arr.reshape((int(rows),int(cols))))
        
        # Parse the header using a loop and
        # the built-in linecache module
        hdr = [getline(source2, i) for i in range(1,7)]
        values = [float(h.split(" ")[-1].strip()) \
         for h in hdr]
        cols,rows,lx,ly,cell,nd = values
        
        if np.all(np.isnan(arr)):
        
            continue        

        if "_vt_" not in file:
            # compute the range of values to plot
            min_arr = np.nanmin(arr)
            max_arr = np.nanmax(arr)
        
            min_int = np.minimum(np.floor(np.log10(max_arr))-4,np.maximum(-5,np.floor(np.log10(min_arr))))
            max_int = np.minimum(1,np.floor(np.log10(max_arr)))
        
        if "_T_" in file:

            # print('temperature')
            levels = np.linspace(min_arr_T,max_arr_T,10)
            label_str = 'Temperature (K)'

        elif "_dep" in file:

            # print('deposit')
            levels = 10.0**np.linspace(min_int_dep,max_int_dep,2*int(max_int_dep-min_int_dep)+1)
            label_str = 'Deposit thickness (m)'

        elif "_ers" in file:

            # print('erosion')
            levels = 10.0**np.linspace(min_int_ers,max_int_ers,2*int(max_int_ers-min_int_ers)+1)
            label_str = 'Erosion thickness (m)'

        elif "_vt_" in file:

            # print('dynamic pressure')
            levels = dyn_pressure_unique
            label_str = 'Dynamic pressure (Pa)'
            lev = 0
            for dyn_lev in dyn_pressure_unique:
                arr[arr==lev]=dyn_lev
                lev +=1

            levels = np.append(levels,np.max(levels)*2)

        else:

            # print('thickness')
            levels = 10.0**np.linspace(min_int_thk,max_int_thk,2*int(max_int_thk-min_int_thk)+1)
            label_str = 'Flow thickness (m)'

        B = np.ma.array(arr, mask=np.isnan(arr))
        
        # ravel() converts C to a 1d-array
        # p1.set_array(B.ravel())
        p1.set_array(B)

        
        cmap = plt.get_cmap('terrain_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        p1.set_norm(norm=norm)

        clb.remove()


        if "_vt_" in file:
            ticks  = 0.5*(levels[:-1]+levels[1:])
            clb = plt.colorbar(format=ticker.FuncFormatter(fmt),ticks=ticks)
            clb.set_ticklabels(np.char.mod('%d', levels[:-1]))
            
        else:
        
            clb = plt.colorbar(format=ticker.FuncFormatter(fmt))

        clb.set_label(label_str, labelpad=-40, y=1.05, rotation=0)

        txt_box.remove()
        
        txt_box = plt.text(0.1, 0.9,t_string,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes,bbox={'facecolor': 'white', 'alpha': 0.7, 'pad': 5})            
        
        plt.axis('scaled')
        plt.xlim([plt_xmin,plt_xmax])
        plt.ylim([plt_ymin,plt_ymax])
        plt.tight_layout()
        plt.savefig(source2.replace('asc','jpg'),dpi=200)

        
    counter+=1
        



