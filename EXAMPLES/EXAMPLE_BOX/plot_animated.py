#!/usr/bin/env python
"""
Animate the 1D output of IMEX_SfloW2D

"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

runname = sys.argv[1]

millisec = sys.argv[2]

bakfile = runname + '.bak'

with open(bakfile) as f:
    for line in f:
        if "DT_OUTPUT" in line:
            line1 = line.replace(',', '')
            dt = float(line1.split()[-1])
            print('dt=' + str(dt))

filelist = sorted(glob.glob(runname + '_*[0-9].p_2d'))

filename = runname + '_{0:04}'.format(0) + '.p_2d'
data = np.loadtxt(filename, skiprows=0)

x_cent = data[:, 0]
y_cent = data[:, 1]
h = data[:, 2]
u = data[:, 3]
v = data[:, 4]
B_cent = data[:, 5]
hB = data[:, 6]
xs = data[:, 7]

x_unique = np.unique(x_cent)

n_cent = len(x_cent)
n_unique = len(x_unique)

# create a figure for the plot
fig, ax = plt.subplots()

plt.xlim([np.amin(x_unique), np.amax(x_unique)])
plt.ylim([-1, 1000])
# plt.gca().set_aspect('equal', adjustable='box')

line1, = ax.plot(x_unique, B_cent)
line2, = ax.plot(x_unique, hB)

time_template = 'time = %.2fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

ax.grid()


def animate(i):

    filename = runname + '_{0:04}'.format(i) + '.p_2d'
    print(filename)
    data = np.loadtxt(filename, skiprows=0)

    x_cent = data[:, 0]
    B_cent = data[:, 5]
    hB = data[:, 6]

    x_unique = np.unique(x_cent)

    line1.set_data(x_unique, B_cent)
    line2.set_data(x_unique, hB)
    time_text.set_text(time_template % (i * dt))

    return line2, time_text


ani = animation.FuncAnimation(fig,
                              animate,
                              np.arange(0, len(filelist)),
                              interval=millisec)

ani.save(runname + '.mp4', fps=5)
plt.show()
