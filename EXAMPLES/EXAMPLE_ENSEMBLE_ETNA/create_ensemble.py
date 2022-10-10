import pandas as pd
from scipy.stats import norm
import numpy as np
import random
from linecache import getline


##########################################
def uniform_discrete_2D_sampling(seq_1, seq_2, n_samples):

    # Sample from arrays

    list_indx = random.choices(range(len(seq_1)), k=n_samples)

    list_1 = [seq_1[i] for i in list_indx]
    list_2 = [seq_2[i] for i in list_indx]

    return list_1, list_2


##########################################
def log_continuous_1D_sampling(val_min, val_max, n_samples):

    # Sample from interval with log distribution

    log_val = np.random.uniform(np.log(val_min), np.log(val_max), n_samples)
    val = np.exp(log_val)

    return val


##########################################
def linear_continuous_1D_sampling(val_min, val_max, n_samples):

    # Sample from interval with uniform distribution

    val = np.random.uniform(val_min, val_max, n_samples)

    return val


##########################################

def associate_sector_to_coords(sect_asc_file,coor_x,coor_y):

    # Read .asc file into a numpy array
    print('Reading .asc file: ' + sect_asc_file)

    hdr = [getline(sect_asc_file, i) for i in range(1, 7)]
    values = [float(h.split(" ")[-1].strip())
              for h in hdr]

    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    val = pd.read_table(sect_asc_file,
                        delim_whitespace=True,
                        header=None,
                        skiprows=6,
                        dtype='unicode').astype(float).values
    val = np.flipud(val)  
    
    x_min = lx
    x_max = lx + cols * cell

    y_min = lx
    y_max = ly + rows * cell
          
    sector_list = []
    
    for xp,yp in zip(coor_x,coor_y):
    
        # get the index of the cell from x and y
        ix = int(np.floor((xp - lx) / cell))
        ix = max(0,min(ix,cols-1))
        iy = int(np.floor((yp - ly) / cell))
        iy = max(0,min(iy,rows-1))
        
        sector_list.append(int(val[iy,ix]))
      
    return sector_list
    
    
##########################################
def nonuniform_continuous_grid_sampling(asc_file, n_samples, discrete_flag):

    # Sample from non-uniform distribution on uniform 1D/2D grid

    # Read .asc file into a numpy array
    print('Reading .asc file: ' + asc_file)

    hdr = [getline(asc_file, i) for i in range(1, 7)]
    values = [float(h.split(" ")[-1].strip())
              for h in hdr]

    cols, rows, lx, ly, cell, nd = values
    cols = int(cols)
    rows = int(rows)

    val = pd.read_table(asc_file,
                        delim_whitespace=True,
                        header=None,
                        skiprows=6,
                        dtype='unicode').astype(float).values
    val = np.flipud(val)

    # normalize the values
    val = val / np.sum(val)

    # print the normalized values
    # print(val)

    # define the ranges for the variables
    val_max = np.amax(val)

    x_min = lx
    x_max = lx + cols * cell

    y_min = lx
    y_max = ly + rows * cell

    # initialize the list for the two sample arrays
    x_list = []
    y_list = []

    for i in range(n_samples):
        # loop to build the sample

        while True:
            # iterate until a point (x,y) is valid

            # generate a random point in the 3D space (x,y,val)
            x_test = np.random.uniform(x_min, x_max)
            y_test = np.random.uniform(y_min, y_max)
            val_test = np.random.uniform(0, val_max)

            # get the index of the cell from x and y
            ix = int(np.floor((x_test - lx) / cell))
            iy = int(np.floor((y_test - ly) / cell))

            if (val_test <= val[iy, ix]):

                if discrete_flag:

                    # if we want the point in the pixel center
                    x_list.append(x_min + cell * (0.5 * ix))
                    y_list.append(y_min + cell * (0.5 * iy))

                else:

                    # if (x,y,val) is below the step-surface defined by
                    # the .asc file accept the point (x,y) and search
                    # for a new one. The probability to be below the
                    # surface is proportional to val[i,j]
                    x_list.append(x_test)
                    y_list.append(y_test)

                break

    # Compute the bi-dimensional histogram of the two data samples x_asc and y_asc
    # The histogram is normalized, and it should converge to the normalized ascii
    # file.
    H, yedges, xedges = np.histogram2d(x_list,
                                       y_list,
                                       bins=(cols, rows),
                                       range=[[x_min, x_max], [y_min, y_max]],
                                       normed=True,
                                       weights=None,
                                       density=None)

    # print the 2D histogram values from the sample
    # print('')
    # print(H.T)

    return x_list, y_list


##########################################
def nonuniform_discrete_sampling(params, weights, n_samples):
    # Sample from list with prescribed weights

    x_min = 0.0
    x_max = len(params)

    w_max = np.amax(weights)

    # initialize the list for the two sample arrays
    param_list = []

    for i in range(n_samples):
        # loop to build the sample

        while True:
            # iterate until a point (x,w) is valid

            # generate a random point in the 3D space (x,w)
            x_test = np.random.uniform(x_min, x_max)
            w_test = np.random.uniform(0, w_max)

            # get the index of the cell from x
            ix = int(x_test - x_min)

            if (w_test <= weights[ix]):

                # if (x,w) is below the step-surface defined by
                # the weights accept the point (x,w) and search
                # for a new one. The probability to be below the
                # surface is proportional to weights[i]
                param_list.append(params[ix])

                break

    return param_list


##########################################

def main():

    list_values = []
    list_names = []

    n_samples = 25

    ###
    """
    seq_x = [500085,500346,499861,500492,500685]
    seq_y = [4177539,4177504,4177572,4177858,4177729]

    list_x , list_y = uniform_discrete_2D_sampling(seq_x,seq_y,n_samples)

    list_values.append(list_x)
    list_names.append('Coor_x')

    list_values.append(list_y)
    list_names.append('Coor_y')
    """

    # range of volumes
    volume_min = 1000000  # m3
    volume_max = 1000000  # m3

    volume = log_continuous_1D_sampling(volume_min, volume_max, n_samples)

    list_values.append(volume)
    list_names.append('Volume')

    ###
    mu_min = 0.1
    mu_max = 0.3

    mu = linear_continuous_1D_sampling(mu_min, mu_max, n_samples)

    list_values.append(mu)
    list_names.append('mu')

    ###
    xi_min = 300.0
    xi_max = 5000.0

    xi = log_continuous_1D_sampling(xi_min, xi_max, n_samples)

    list_values.append(xi)
    list_names.append('xi')

    ###

    """
    asc_file = 'test.asc'
    x_list, y_list = nonuniform_continuous_grid_sampling(asc_file,
                                                         n_samples,
                                                         discrete_flag=False)

    list_values.append(x_list)
    list_names.append('x')

    list_values.append(y_list)
    list_names.append('y')
    """

    ###
    # list of 2D points
    points = [[500360, 4177600], [500650, 4177913]]
    # weight of 2D points
    weights = [2, 1]

    param_list = nonuniform_discrete_sampling(points, weights, n_samples)
    # print(param_list)

    coor_x = [param[0] for param in param_list]
    list_values.append(coor_x)
    list_names.append('Coor_x')

    coor_y = [param[1] for param in param_list]
    list_values.append(coor_y)
    list_names.append('Coor_y')

    ###
    
    """
    sector_list = associate_sector_to_coords('aree.asc',coor_x,coor_y)

    list_values.append(sector_list)
    list_names.append('Sector')
    """
        
    ###


    df = pd.DataFrame(list(map(list, zip(*list_values))), columns=list_names)

    df.to_csv('samples.csv')


if __name__ == '__main__':

    main()
