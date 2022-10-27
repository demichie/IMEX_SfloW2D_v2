import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import os.path


if len(sys.argv)==2: 
 
    prb_file = sys.argv[1]

    os.path.isfile(prb_file) 


else:

    print('Please provide file name (*.prb)\n')
    sys.exit()


data = np.genfromtxt(prb_file, delimiter=',', skip_header = 3)
print(data.shape)


# read the content of the file opened
file = open(prb_file)
content = file.readlines()
  
# read 10th line from the file
# print("variable names line")

vars = content[1].split(',')[:]

# print('vars',vars)

# print("variable units line")

units = content[2].split(',')[:]


for i in range(len(vars)-1):

    # print(i,vars[i+1])

    plt.subplots(figsize=(10, 6))
    plt.plot(data[:,0],data[:,i+1])
    plt.xlim((np.amin(data[:,0]),np.amax(data[:,0])))
    
    plt.xlabel('Time (s)')
    plt.ylabel(vars[i+1]+' '+units[i+1])
    # print((vars[i+1].strip()).split(' ')[0]+'.png')
    fig_name = prb_file.replace('.csv','_'+(vars[i+1].strip().replace('.','_'))+'.png')
    fig_name = fig_name.replace('_.','.')
    plt.savefig(fig_name)

# plt.show()

