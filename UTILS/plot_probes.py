import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path

if len(sys.argv)==2: 
 
    prb_file = sys.argv[1]

    os.path.isfile(prb_file) 


else:

    print('Please provide file name (*_prb_xxxx.csv)\n')
    sys.exit()


# Using readline()
file_prb = open(prb_file, 'r')

line1 = file_prb.readline()
print('line1',line1)
coors = line1.split()
print('coors',coors)
line = file_prb.readline()
vars = line.split(',')
nvars = len(vars)
line = file_prb.readline()
units = line.split(',')

file_prb.close()

data = np.genfromtxt(prb_file, delimiter=',', skip_header = 3)
print(data.shape)

for i in range(nvars-1):

    fig = plt.subplots()

    plt.plot(data[:,0],data[:,i+1])
    plt.xlim([np.amin(data[:,0]),np.amax(data[:,0])])
    plt.xlabel(vars[0].strip()+' '+units[0].strip())
    plt.ylabel(vars[i+1].strip()+' '+units[i+1].strip())
    plt.title('(' + coors[0].strip() + ',' + coors[1].strip() + ')' )
    figfile = prb_file.replace('.csv','_'+str(i+1)+'.pdf')
    plt.savefig(figfile)

# plt.show()

