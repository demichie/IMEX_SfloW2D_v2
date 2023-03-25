import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path

if len(sys.argv) == 2:

    prb_file = sys.argv[1]
    os.path.isfile(prb_file)
    rate = False

elif len(sys.argv) == 3:

    prb_file = sys.argv[1]
    os.path.isfile(prb_file)

    if sys.argv[2] == 'True':

        rate = True

    else:

        rate = False

else:

    print('Please provide file name (*.prb)\n')
    sys.exit()

data = np.genfromtxt(prb_file, delimiter=',', skip_header=3)
# print(data.shape)

# read the content of the file opened
file = open(prb_file)
content = file.readlines()

# read 10th line from the file
# print("variable names line")

vars = content[1].split(',')[:]

# print('vars',vars)

# print("variable units line")

units = content[2].split(',')[:]
print()

for i in range(len(vars) - 1):

    # print(i,vars[i+1])

    plt.rcParams.update({'font.size': 14})
    fig, ax1 = plt.subplots(figsize=(10, 4))
    ax1.plot(data[:, 0], data[:, i + 1])
    ax1.set_xlim((np.amin(data[:, 0]), np.amax(data[:, 0])))

    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel(vars[i + 1].strip() + ' ' + units[i + 1].strip())

    plt.tight_layout()

    fig_name = prb_file.replace(
        '.csv', '_' + (vars[i + 1].strip().replace('.', '_')) + '.pdf')
    fig_name = fig_name.replace('_.', '.')
    plt.savefig(fig_name)

    if rate:

        ax2 = ax1.twinx()
        dx = np.diff(data[:, 0])
        dy = np.diff(data[:, i + 1])
        x_half = 0.5 * (data[:-1, 0] + data[1:, 0])
        dy_dx = dy / dx

        ax2.plot(x_half, dy_dx, 'C1--')
        label = 'Rate ' + units[i + 1].strip().replace(')',
                                                       '') + r'$\cdot s^{-1}$)'
        ax2.set_ylabel(label)
        ax2.tick_params(axis='y', color='C1', labelcolor='C1')

        fig_name = fig_name.replace('.pdf', '_B.pdf')
        plt.tight_layout()
        plt.savefig(fig_name)

# plt.show()
