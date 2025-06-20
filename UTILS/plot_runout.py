import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path

if len(sys.argv) == 2:

    runout_file = sys.argv[1]

    os.path.isfile(runout_file)

else:

    print('Please provide file name (*_runout.csv)\n')
    sys.exit()

data = np.genfromtxt(runout_file, delimiter=',', skip_header=1)
print(data.shape)

time = data[:, 0]
runout = data[:, 1]
area = data[:, 2]

# read the content of the file opened
file = open(runout_file)
content = file.readlines()

# read 10th line from the file
# print("variable names line")

fig, ax1 = plt.subplots()

plt.xlim([np.amin(time), np.amax(time)])

ax2 = ax1.twinx()
ax1.plot(time, runout, 'g-')
ax2.plot(time, area, 'b-')

ax1.set_xlabel('time [s]')
ax1.set_ylabel('runout [m]', color='g')
ax2.set_ylabel('area [m2]', color='b')

fig.savefig(runout_file.replace('csv', 'png'))

idx = int(np.argmax(runout))
print('time', time[idx])
print('runout', runout[idx])

max_runout = runout[idx]

set0 = np.array(np.argwhere(runout > 0.99 * max_runout))
idx = set0[0]
print('99% runout', runout[idx])
print('time for 99% max runout', time[idx])

set0 = np.array(np.argwhere(runout > 0.95 * max_runout))
idx = set0[0]
print('95% runout', runout[idx])
print('time for 95% max runout', time[idx])

# plt.show()
