# defining the libraries
import linecache as lc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import glob
from matplotlib import rc
from matplotlib.colors import LightSource
from matplotlib.colors import BoundaryNorm
import netCDF4
from osgeo import osr


def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print("Folder created at:", folder_path)
    else:
        print("Folder already exists at:", folder_path)


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def fmt2(x, pos):

    return r'${}\%$'.format(x)


def bilin_interp(xin, yin, fin, xout, yout):
    """
    Bilinear interpolation from a regular grid to a list of points

    Parameters:
    - xin, yin, fin: data to be interpolated
    - xout, yout: list of points where to interpolate
    """
    dXin = xin[1] - xin[0]
    dYin = yin[1] - yin[0]

    ix = np.maximum(0, np.ceil((xout - xin[0]) / dXin).astype(int) - 1)
    iy = np.maximum(0, np.ceil((yout - yin[0]) / dYin).astype(int) - 1)

    alfa_x = (xout - xin[ix]) / dXin
    alfa_y = (yout - yin[iy]) / dYin

    fout = alfa_x * alfa_y * fin[iy + 1, ix + 1] + alfa_x * (
        1.0 - alfa_y) * fin[iy, ix + 1] + (1.0 - alfa_x) * alfa_y * fin[
            iy + 1, ix] + (1.0 - alfa_x) * (1.0 - alfa_y) * fin[iy, ix]

    return fout


def save_to_csv(file_path, column1, column2, column3, column4):
    """
    Save data from four columns into a CSV file.

    Parameters:
    - file_path: The path to the CSV file.
    - column1, column2, column3, column4: Lists containing data for each
      column.
    """

    data = {
        'Water Depth [m]': column1,
        '5%ile': column2,
        '50%ile': column3,
        '95%ile': column4
    }
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    df.to_csv(file_path, index=False)
    print(f'Data saved to {file_path}')


def read_asc_file(file_path):
    """
    Read data from a ESRI ascii file.

    Parameters:
    - file_path: The path to the .asc file.
    """

    try:
        # Parse the header using a loop and
        # the built-in linecache module
        hdr = [lc.getline(file_path, i) for i in range(1, 7)]
        values = [float(h.split(" ")[-1].strip()) for h in hdr]
        ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value = values
        ncols = int(ncols)
        nrows = int(nrows)

        column_names = [i for i in range(0, ncols)]

        data = pd.read_csv(file_path,
                           delim_whitespace=True,
                           skiprows=6,
                           header=None,
                           nrows=nrows,
                           names=column_names).to_numpy()

        return data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value

    except Exception as e:
        print(f"Error reading ASC file: {e}")
        return None


rc("pdf", fonttype=42)


# Print iterations progress
def printProgressBar(iteration,
                     total,
                     prefix='',
                     suffix='',
                     decimals=1,
                     bar_length=20):
    """
    Call in a loop to create terminal progress bar

    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' %
                     (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def main(argv):

    EPSG = 32633

    h_thresholds = [0.1, 0.5, 1.0, 2.0, 5.0]

    current_dir = os.getcwd()

    # Use glob to find folders containing the string "ensemble"
    ensemble_folders = glob.glob(os.path.join(current_dir, '*ensemble*'))

    # Filter out only directories from the results
    ensemble_folders = [
        folder for folder in ensemble_folders if os.path.isdir(folder)
    ]

    folder = ensemble_folders[0]
    print(folder)
    fl = glob.glob(folder + '/dem_esri.asc')
    print(fl[0])

    h, topo_nx, topo_ny, topo_xll, topo_yll, topo_cellsize, topo_nodata_value = read_asc_file(
        fl[0])

    xmin = topo_xll
    xmax = topo_xll + topo_nx * topo_cellsize

    ymin = topo_yll
    ymax = topo_yll + topo_ny * topo_cellsize

    print(xmin, xmax, ymin, ymax)

    folder = ensemble_folders[0]
    fl = glob.glob(folder + '/*_max.asc')
    print(fl[0])

    asc0, nx, ny, xref, yref, cellsize, nodata_value = read_asc_file(fl[0])

    extent1 = [xref, xref + nx * cellsize, yref, yref + ny * cellsize]

    X_grid = np.zeros((ny, nx))

    shape_hMax = X_grid.shape + tuple([len(ensemble_folders)])

    hMax = np.zeros(shape_hMax)
    # print(wd_np.shape)

    print('Reading maps')

    for i, folder in enumerate(ensemble_folders):

        printProgressBar(i, len(ensemble_folders) - 1)

        # print(folder)
        fl = glob.glob(folder + '/*_max.asc')
        # print(fl)

        asc, nx, ny, xref, yref, dxy, nodata_value = read_asc_file(fl[0])

        asc[asc == nodata_value] = 0.0

        hMax[..., i] = asc

    hMax_min = np.min(hMax, axis=-1)
    hMax_max = np.max(hMax, axis=-1)

    # 2d array for pixels where there are are different values of wd
    check_array = hMax_min != hMax_max

    h_thresholds = [0.1, 0.5, 1.0, 5.0]
    n_thr = len(h_thresholds)

    probs = np.zeros((n_thr, ny, nx))

    print('Processing maps')

    hMax_ij = np.zeros(len(ensemble_folders))

    for j in range(nx):

        printProgressBar(j, nx - 1)

        for i in range(ny):

            if check_array[i, j]:

                hMax_ij = hMax[i, j, :]

                for idx, threshold in enumerate(h_thresholds):

                    above_threshold = hMax_ij > threshold
                    count = np.sum(above_threshold)
                    probs[idx, i, j] = count * 100.0 / len(ensemble_folders)

            else:

                probs[:, i, j] = 0

    output_path = "./exceedance_maps"
    create_folder_if_not_exists(output_path)
    
    output_folder = glob.glob(os.path.join(current_dir,'exceedance_maps'))[0]
    print(output_folder)

    header = "ncols     %s\n" % nx
    header += "nrows    %s\n" % ny
    header += "xllcorner " + str(xref) + "\n"
    header += "yllcorner " + str(yref) + "\n"
    header += "cellsize " + str(dxy) + "\n"
    header += "NODATA_value " + str(nodata_value)

    x_stag = np.linspace(xref, xref + nx * dxy, int(nx + 1))
    y_stag = np.linspace(yref, yref + ny * dxy, int(ny + 1))

    X, Y = np.meshgrid(x_stag, y_stag)

    ncfilename = output_folder + '/exceedance_probability.nc'

    print(ncfilename)

    ncfile = netCDF4.Dataset(ncfilename, mode='w', format='NETCDF4_CLASSIC')
    ncfile.createDimension('x', nx)
    ncfile.createDimension('y', ny)
    # unlimited axis (can be appended to).
    ncfile.createDimension('h_thr', None)

    ncfile.title = 'hMax percentiles output'

    ncfile.Conventions = "CF-1.0"
    ncfile.subtitle = "My model data subtitle"
    ncfile.anything = "write anything"

    x = ncfile.createVariable('x', np.float64, ('x', ))
    x.long_name = 'Easting [UTM]'
    x.units = 'meters'
    x[:] = x_stag[:-1] + 0.5 * dxy

    y = ncfile.createVariable('y', np.float64, ('y', ))
    y.long_name = 'Northing [UTM]'
    y.units = 'meters'
    y[:] = y_stag[:-1] + 0.5 * dxy

    h_thr = ncfile.createVariable('h_thr', np.float64, ('h_thr', ))
    h_thr.long_name = 'Threshold'
    h_thr.units = 'meters'

    prob_hMax = ncfile.createVariable('prob_hMax',
                                      np.float64, ('h_thr', 'y', 'x'),
                                      zlib=True,
                                      fill_value=nodata_value)
    prob_hMax.standard_name = 'probability'
    prob_hMax.units = 'percentage'

    # https://gis.stackexchange.com/questions/373193/
    #         create-utm-netcdf-to-be-opened-in-arcgis-pro
    proj = osr.SpatialReference()
    proj.ImportFromEPSG(EPSG)

    prob_hMax.esri_pe_string = proj.ExportToWkt()

    for i, threshold in enumerate(h_thresholds):

        h_thr[i] = threshold

        asc_file = output_folder + '/hazard_map_' + str(threshold).zfill(2) + '.asc'
        png_file = output_folder + '/hazard_map_' + str(threshold).zfill(2) + '.png'

        print(asc_file)

        asc = np.squeeze(probs[i, :, :])
        prob_hMax[i, :, :] = np.flipud(asc)

        asc[asc == 0] = nodata_value

        np.savetxt(asc_file, asc, header=header, fmt='%1.5f', comments='')

        fig, ax = plt.subplots()
        fig.set_size_inches(2.0 + 6.0 / ny * nx, 6.0)

        ls = LightSource(azdeg=315, altdeg=45)

        plt.imshow(ls.hillshade(np.flipud(h),
                                vert_exag=1.0,
                                dx=cellsize,
                                dy=cellsize),
                   cmap='gray',
                   extent=extent1)

        im_ratio = ny / nx
        cmap = plt.get_cmap('viridis')
        # levels = np.linspace(0, 10, 11)
        levels = [0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        label_str = 'Exceedence probability [%]'

        norm = BoundaryNorm(levels, ncolors=cmap.N)
        asc[asc == nodata_value] = np.nan
        p1 = plt.pcolormesh(X,
                            Y,
                            np.flipud(asc),
                            cmap=cmap,
                            norm=norm,
                            alpha=0.65)

        clb = plt.colorbar(
            p1,
            fraction=0.046 * im_ratio,
            pad=0.04,
            # format=ticker.FuncFormatter(fmt),
            ax=ax)

        clb.set_label(label_str)
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])
        plt.tight_layout()
        plt.savefig(png_file, dpi=600)
        # plt.show()
        plt.close(fig)

    print(ncfile)
    ncfile.close()
    print('Dataset is closed!')


if __name__ == "__main__":
    main(sys.argv[1:])
