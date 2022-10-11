import numpy as np
import pandas as pd
from linecache import getline
import os
import csv


def calculate_map(asc_file, row_list, group_index):

    n_rows = len(row_list)
    FLAG = "FALSE"  # necessary for initialization of array sum
    c = 0  # variable counting not existing file when the simulation crashed

    for i in row_list:

        # changing path to the ensemble.xxxxx dir containing .asc file

        working_dir = "ensemble." + "{:05d}".format(i)
        working_dir = os.path.join(os.getcwd(), working_dir)
        main_dir = os.getcwd()
        os.chdir(working_dir)

        if os.path.isfile(asc_file):

            # reading .asc file

            print('Reading .asc file: ' + asc_file)

            hdr = [getline(asc_file, i) for i in range(1, 7)]
            values = [float(h.split(" ")[-1].strip()) for h in hdr]

            cols, rows, lx, ly, cell, nd = values
            cols = int(cols)
            rows = int(rows)

            array_inv = pd.read_table(asc_file,
                                      delim_whitespace=True,
                                      header=None,
                                      skiprows=6,
                                      dtype='unicode').astype(float).values
            array_inv = np.flipud(array_inv)

            ######

            if FLAG == "FALSE":
                # array initialization which contain the sum of array invasion
                array_inv_tot = np.zeros_like(array_inv)
                FLAG = "TRUE"  # initialization ok

            # sum of asc file with 0 and 1 values: array_inv>0 return True (or 1) if values are positive, otherwise return FALSE (zero)
            array_inv_tot += (array_inv > 0)

        else:

            print('No found .asc file: ' + asc_file)
            print(working_dir)
            c += 1  # counting not existing file

        os.chdir(main_dir)  # coming back to the main directory

    array_inv_tot = array_inv_tot / (n_rows - c)  # array_inv_tot normalization

    # dir containing new asc files
    result_dir = "result." + "{:03d}".format(group_index)
    result_dir = os.path.join(os.getcwd(), result_dir)

    isExist = os.path.exists(result_dir)

    if not isExist:
        os.mkdir(result_dir)  # creates the dir if it not exist

    os.chdir(result_dir)

    header = "ncols     %s\n" % cols
    header += "nrows    %s\n" % rows
    header += "xllcorner " + str(lx) + "\n"
    header += "yllcorner " + str(ly) + "\n"
    header += "cellsize " + str(cell) + "\n"
    header += "NODATA_value " + str(nd)

    np.savetxt(asc_file,
               np.flipud(array_inv_tot),
               header=header,
               fmt='%1.5f',
               comments='')

    os.chdir(main_dir)


#####  main  #####
def main():

    current_dir = os.getcwd()
    os.chdir('./templatedir/')

    with open('IMEX_SfloW2D.template') as fp:

        thick_levels = []
        dyn_pres_levels = []

        for cnt, line in enumerate(fp):

            if "RUN_NAME" in line:

                run_name = line.replace('RUN_NAME', '')
                run_name = run_name.replace('=', '')
                run_name = run_name.replace('"', '')
                run_name = run_name.replace(',', '')
                run_name = run_name.replace(' ', '')
                run_name = run_name.rstrip()

            if "THICKNESS_LEVELS0" in line:

                thick_levels = line.replace('THICKNESS_LEVELS0', '')
                thick_levels = thick_levels.replace('=', '')
                thick_levels = thick_levels.replace('"', '')
                thick_levels = thick_levels.split(',')
                thick_levels = [
                    string.replace(',', '') for string in thick_levels
                ]
                thick_levels = [
                    string.replace(' ', '') for string in thick_levels
                ]
                thick_levels = [string.rstrip() for string in thick_levels]
                thick_levels = list(filter(None, thick_levels))
                # print(thick_levels)
                # print(len(thick_levels))

            if "DYN_PRES_LEVELS0" in line:

                dyn_pres_levels = line.replace('DYN_PRES_LEVELS0', '')
                dyn_pres_levels = dyn_pres_levels.replace('=', '')
                dyn_pres_levels = dyn_pres_levels.replace('"', '')
                dyn_pres_levels = dyn_pres_levels.split(',')
                dyn_pres_levels = [
                    string.replace(',', '') for string in dyn_pres_levels
                ]
                dyn_pres_levels = [
                    string.replace(' ', '') for string in dyn_pres_levels
                ]
                dyn_pres_levels = [
                    string.rstrip() for string in dyn_pres_levels
                ]
                dyn_pres_levels = list(filter(None, dyn_pres_levels))
                # print(dyn_pres_levels)
                # print(len(dyn_pres_levels))

    n_VT = len(thick_levels) * len(dyn_pres_levels)
    print('n_VT', n_VT)

    os.chdir(current_dir)

    # read the csv generated by the ensemble by sampling vents position, volumes, parameters
    df = pd.read_csv('samples.csv')

    # Single group: uncomments if grouping is not used
    n_groups = 1
    row_list = range(len(df.index))

    print(row_list)

    fname = run_name + "_max.asc"

    calculate_map(fname, row_list, 0)

    for i in range(n_VT):
        VTfile = run_name + "_VT" + "_{:04d}".format(i + 1)
        VTfile = VTfile + ".asc"
        # calculate map invasion for each combination of thickness-dynamic pressure threshold
        calculate_map(VTfile, row_list, 0)


    # group by coordinates (both x and y)
    group = df.groupby(["Coor_x", "Coor_y"])
    n_groups = len(group.groups)
    group_list = list(group.groups.values())
    """

    # create the bins (intervals) for mu
    nbins_mu = 10 
    mu_bins = np.linspace(0.1, 0.3, nbins_mu)

    # create the bins (intervals) for xi
    nbins_xi = 10
    xi_bins = np.linspace(300.0,5000.0, nbins_xi)

    # group the dataframe elements in sub-domains by using the two partitions
    groupBymu = df.groupby(pd.cut(df['mu'], mu_bins))
    groupByxi = df.groupby(pd.cut(df['xi'], xi_bins))
    groupsByxiAndmu = df.groupby([pd.cut(df['mu'], mu_bins), pd.cut(df['xi'], xi_bins)])

    """

    fname = run_name + "_max.asc"

    for i_group in range(n_groups):

        if n_groups > 1:

            group_file = 'samples' + "_{:03d}".format(i_group + 1) + '.csv'
            row_list = group_list[i_group].tolist()
            df.iloc[row_list].to_csv(group_file, index=False)
        
        calculate_map(fname, row_list, i_group+1)

        for i in range(n_VT):
            VTfile = run_name + "_VT" + "_{:04d}".format(i + 1)
            VTfile = VTfile + ".asc"
            # calculate map invasion for each combination of thickness-dynamic pressure threshold
            calculate_map(VTfile, row_list, i_group+1)

        

if __name__ == '__main__':

    main()
