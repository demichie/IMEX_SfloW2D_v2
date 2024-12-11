import matplotlib.pyplot as plt
import numpy as np
from linecache import getline
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as ticker
import glob
import os
from matplotlib.colors import LightSource

topography_file = '../DEM/Etna2014_crop.asc'


def interp2Dgrids(xin, yin, Zin, Xout, Yout):
    """
    Interpolation from a regular grid to a second regular grid

    @params:
        xin      - Required : original grid X values (1D Dble)
        yin      - Required : original grid Y values (1D Dble)
        Zin      - Required : original grid Z values (2D Dble)
        xout     - Required : new grid X values (2D Dble)
        yout     - Required : new grid Y values (2D Dble)
    """
    xinMin = np.min(xin)

    yinMin = np.min(yin)

    cellin = xin[1] - xin[0]

    xout = Xout[0, :]

    yout = Yout[:, 0]

    # Search for the cell
    xi = (xout - xinMin) / cellin
    yi = (yout - yinMin) / cellin

    # Indexes of the lower-left corner of the cell
    ix = np.maximum(0, np.minimum(xin.shape[0] - 2, np.floor(xi).astype(int)))
    iy = np.maximum(0, np.minimum(yin.shape[0] - 2, np.floor(yi).astype(int)))

    # Indexes of the top-right corner of the cell
    ix1 = ix + 1
    iy1 = iy + 1

    # Relative coordinates of the center of the parent lobe in the cell
    xi_fract = np.maximum(0, np.minimum(1, (xi - ix).reshape(1,
                                                             Xout.shape[1])))
    yi_fract = np.maximum(0, np.minimum(1, (yi - iy).reshape(Yout.shape[0],
                                                             1)))

    xi_out_yi = np.outer(yi_fract, xi_fract)

    Zout = xi_out_yi * Zin[np.ix_(iy1, ix1)] + \
        (xi_fract - xi_out_yi) * Zin[np.ix_(iy, ix1)] + \
        (yi_fract - xi_out_yi) * Zin[np.ix_(iy1, ix)] + \
        (1.0 - xi_fract - yi_fract + xi_out_yi) * Zin[np.ix_(iy, ix)]

    return Zout


def regrid(xin, yin, fin, xl, xr, yl, yr):

    nXin = xin.shape[0] - 1
    nYin = yin.shape[0] - 1

    dXin = xin[1] - xin[0]
    dYin = yin[1] - yin[0]

    ix1 = int(np.maximum(0, np.ceil((xl - xin[0]) / dXin).astype(int) - 1))
    ix2 = int(np.minimum(nXin, np.ceil((xr - xin[0]) / dXin).astype(int)))

    iy1 = int(np.maximum(0, np.ceil((yl - yin[0]) / dYin).astype(int) - 1))
    iy2 = int(np.minimum(nYin, np.ceil((yr - yin[0]) / dYin).astype(int)))

    print('ix1,ix2,iy1,iy2', ix1, ix2, iy1, iy2)

    fout = 0.0

    for ix in range(ix1, ix2):

        alfa_x = (np.minimum(xr, xin[ix + 1]) -
                  np.maximum(xl, xin[ix])) / (xr - xl)

        for iy in range(iy1, iy2):

            alfa_y = (np.minimum(yr, yin[iy + 1]) -
                      np.maximum(yl, yin[iy])) / (yr - yl)

            fout = fout + alfa_x * alfa_y * fin[iy, ix]

    return fout


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


current_dir = os.getcwd()
os.chdir('./templatedir/')

with open('IMEX_SfloW2D.inp') as fp:

    for cnt, line in enumerate(fp):

        if "RUN_NAME" in line:

            run_name = line.replace('RUN_NAME', '')
            run_name = run_name.replace('=', '')
            run_name = run_name.replace('"', '')
            run_name = run_name.replace(',', '')
            run_name = run_name.replace(' ', '')
            run_name = run_name.rstrip()

        if "T_START" in line:
            t_start_str = line.replace('T_START', '')
            t_start_str = t_start_str.replace('=', '')
            t_start_str = t_start_str.replace(',', '')
            t_start = float(t_start_str)
            # print("t_start =",t_start)

        if "T_END" in line:
            t_end_str = line.replace('T_END', '')
            t_end_str = t_end_str.replace('=', '')
            t_end_str = t_end_str.replace(',', '')
            print("t_end =", t_end_str)
            t_end = float(t_end_str)
            print("t_end =", t_end)

        if "DT_OUTPUT" in line:
            dt_output_str = line.replace('DT_OUTPUT', '')
            dt_output_str = dt_output_str.replace('=', '')
            dt_output_str = dt_output_str.replace(',', '')
            dt_output = float(dt_output_str)
            # print("dt_output =",dt_output)

        if "X0=" in line:
            x0_str = line.replace('X0', '')
            x0_str = x0_str.replace('=', '')
            x0_str = x0_str.replace(',', '')
            x0 = float(x0_str)
            print("x0 =", x0)

        if "Y0=" in line:
            y0_str = line.replace('Y0', '')
            y0_str = y0_str.replace('=', '')
            y0_str = y0_str.replace(',', '')
            y0 = float(y0_str)
            print("y0 =", y0)

        if "COMP_CELLS_X" in line:
            comp_cells_x_str = line.replace('COMP_CELLS_X', '')
            comp_cells_x_str = comp_cells_x_str.replace('=', '')
            comp_cells_x_str = comp_cells_x_str.replace(',', '')
            comp_cells_x = float(comp_cells_x_str)
            print("comp_cells_x =", comp_cells_x)

        if "COMP_CELLS_Y" in line:
            comp_cells_y_str = line.replace('COMP_CELLS_Y', '')
            comp_cells_y_str = comp_cells_y_str.replace('=', '')
            comp_cells_y_str = comp_cells_y_str.replace(',', '')
            comp_cells_y = float(comp_cells_y_str)
            print("comp_cells_y =", comp_cells_y)

        if "CELL_SIZE" in line:
            cell_size_str = line.replace('CELL_SIZE', '')
            cell_size_str = cell_size_str.replace('=', '')
            cell_size_str = cell_size_str.replace(',', '')
            cell_size_str = cell_size_str.replace('D', '')
            cell_size = float(cell_size_str)
            print("cell_size =", cell_size)

print(run_name)
print(len(run_name))

plt_xmin = x0
plt_xmax = x0 + comp_cells_x * cell_size
plt_ymin = y0
plt_ymax = y0 + comp_cells_y * cell_size
print('Topography file ', topography_file)
# Parse the header using a loop and
# the built-in linecache module
hdr = [getline(topography_file, i) for i in range(1, 7)]
values = [float(h.split(" ")[-1].strip()) for h in hdr]
cols, rows, lx, ly, cell1, nd = values
xres = cell1
yres = cell1 * -1

# Load the dem into a numpy array
h = np.loadtxt(topography_file, skiprows=7)
h = np.flipud(h)
# make these smaller to increase the resolution
dx, dy = 0.05, 0.05

x = np.arange(-3.0, 3.0, dx)
y = np.arange(-3.0, 3.0, dy)
X, Y = np.meshgrid(x, y)

xin = np.linspace(lx, lx + (cols) * cell1, int(cols + 1))
yin = np.linspace(ly, ly + (rows) * cell1, int(rows + 1))

# when layering multiple images, the images need to have the same
# extent.  This does not mean they need to have the same shape, but
# they both need to render to the same coordinate system determined by
# xmin, xmax, ymin, ymax.  Note if you use different interpolations
# for the images their apparent extent could be different due to
# interpolation edge effects

extent1 = [lx, lx + cols * cell1, ly, ly + rows * cell1]
"""
fig, ax = plt.subplots(1,1,figsize=(6.0/rows*cols, 6.0))

ls = LightSource(azdeg=315, altdeg=45)


im1 = plt.imshow(ls.hillshade(np.flipud(h),vert_exag=1.0,
                     dx=cell,dy=cell),cmap='gray',extent=extent1)

plt.xlim([plt_xmin, plt_xmax])
plt.ylim([plt_ymin, plt_ymax])

ax.set_aspect('equal', 'box')


# plt.show()
"""

os.chdir(current_dir)

min_arr = 0.0
max_arr = 1.0

for dir in glob.glob("result*"):

    os.chdir(dir)

    for file in glob.glob("*.asc"):

        t_output = t_end

        t_string = 't={t_output:.1f}s'.format(t_output=t_output)
        # print('t_string',t_string,t_output)

        source2 = file
        print(source2)

        # Parse the header using a loop and
        # the built-in linecache module
        hdr = [getline(source2, i) for i in range(1, 7)]
        values = [float(h.split(" ")[-1].strip()) for h in hdr]
        cols, rows, lx, ly, cell, nd = values
        xres = cell
        yres = cell * -1

        # Load the dem into a numpy array
        arr = np.loadtxt(source2, skiprows=6)
        arr[arr == nd] = np.nan
        arr = np.flipud(arr)

        if np.all(np.isnan(arr)):

            continue

        # when layering multiple images, the images need to have the same
        # extent.  This does not mean they need to have the same shape, but
        # they both need to render to the same coordinate system determined by
        # xmin, xmax, ymin, ymax.  Note if you use different interpolations
        # for the images their apparent extent could be different due to
        # interpolation edge effects

        x = np.linspace(lx, lx + (cols) * cell, int(cols + 1))
        y = np.linspace(ly, ly + (rows) * cell, int(rows + 1))

        X, Y = np.meshgrid(x, y)

        extent2 = lx, lx + cols * cell, ly, ly + rows * cell

        levels = np.linspace(min_arr, max_arr, 11)
        label_str = 'Probabilty [0;1]'

        fig, ax = plt.subplots()
        fig.set_size_inches(1.5 + 6.0 / rows * cols, 6.0)

        ls = LightSource(azdeg=45, altdeg=45)

        h_new = interp2Dgrids(xin, yin, h, X, Y)

        im1 = plt.imshow(ls.hillshade(np.flipud(h_new),
                                      vert_exag=2.0,
                                      dx=cell,
                                      dy=cell),
                         cmap='gray',
                         interpolation='nearest',
                         extent=extent1)

        plt.xlim([plt_xmin, plt_xmax])
        plt.ylim([plt_ymin, plt_ymax])

        ax.set_aspect('equal', 'box')

        im_ratio = rows / cols

        cmap = plt.get_cmap('terrain_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        p1 = plt.imshow(arr,
                        cmap=cmap,
                        interpolation='nearest',
                        extent=extent1,
                        alpha=0.65,
                        norm=norm)

        clb = plt.colorbar(p1,
                           fraction=0.046 * im_ratio,
                           pad=0.04,
                           format=ticker.FuncFormatter(fmt),
                           ax=ax)
        clb.set_ticks(np.linspace(0.05, 0.95, num=10))
        ticks = [
            "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6",
            "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"
        ]
        clb.set_ticklabels(ticks)
        clb.set_label(label_str, labelpad=-40, y=1.05, rotation=0)

        txt_box = plt.text(0.1,
                           0.9,
                           t_string,
                           horizontalalignment='center',
                           verticalalignment='center',
                           transform=ax.transAxes,
                           bbox={
                               'facecolor': 'white',
                               'alpha': 0.7,
                               'pad': 5
                           })

        plt.xlim([plt_xmin, plt_xmax])
        plt.ylim([plt_ymin, plt_ymax])
        plt.tight_layout()
        plt.savefig(source2.replace('asc', 'png'), dpi=200)
        plt.close(fig)

    os.chdir(current_dir)
