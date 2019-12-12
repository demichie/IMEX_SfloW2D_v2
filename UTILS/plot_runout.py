#!/usr/bin/env python
"""
% This function 

"""
import numpy as np                      
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import sys
import os.path

print 

if len(sys.argv)==2: 
 
    filename = sys.argv[1]

    os.path.isfile(filename) 


else:

    print('Please provide file name (*_runout.txt)\n')
    sys.exit()

time =[]
area =[]
runout = []

with open(filename) as f:
    lines=f.readlines()
    for line in lines:
        strings = line.split()
        time.append(np.float(strings[3]))
        area.append(np.float(strings[11]))
        runout.append(np.float(strings[7]))
# create a figure for the plot

fig, ax1 = plt.subplots()

plt.xlim([np.amin(time),np.amax(time)])

ax2 = ax1.twinx()
ax1.plot(time, runout, 'g-')
ax2.plot(time, area, 'b-')

ax1.set_xlabel('time [s]')
ax1.set_ylabel('runout [m]', color='g')
ax2.set_ylabel('area [m2]', color='b')

plt.show()

#plt.show()





