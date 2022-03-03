import numpy as np
import matplotlib.pyplot as plt

prb_file = 'exampleKRAKATAU_prb_0005.csv'

data = np.genfromtxt(prb_file, delimiter=',', skip_header = 3)
print(data.shape)


# read the content of the file opened
file = open(prb_file)
content = file.readlines()
  
# read 10th line from the file
# print("variable names line")

vars = content[1].split(',')[:-1]

# print("variable units line")

units = content[2].split(',')[:-1]


for i in range(len(vars)-1):

    # print(i,vars[i+1])

    plt.subplots(figsize=(10, 6))
    plt.plot(data[:,0],data[:,i+1])
    
    plt.xlabel('Time (s)')
    plt.ylabel(vars[i+1]+' '+units[i+1])
    # print((vars[i+1].strip()).split(' ')[0]+'.png')
    fig_name = prb_file.replace('.csv','_'+(vars[i+1].strip())+'.png')
    plt.savefig(fig_name)

# plt.show()

