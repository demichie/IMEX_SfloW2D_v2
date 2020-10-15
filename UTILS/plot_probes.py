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

    print('Please provide file name (*.prb)\n')
    sys.exit()

time =[]
thickness = []

with open(filename) as f:
    lines=f.readlines()
    for line in lines:
        strings = line.split()
        time.append(np.float(strings[0]))
        thickness.append(np.float(strings[1]))
# create a figure for the plot

fig, ax1 = plt.subplots()

plt.xlim([np.amin(time[1:]),np.amax(time[1:])])

ax1.plot(time[1:], thickness[1:], 'b-')

ax1.set_xlabel('time [s]')
ax1.set_ylabel('thickness [m]', color='b')

plt.show()

#plt.show()





