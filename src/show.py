from matplotlib import pyplot as plt
import numpy as np

def read_data(inp):

	retval = list()
	with open(inp, 'r') as fin:
		for line in fin:
			retval.append(float(line))
	
	return np.array(retval)



data = read_data('log')
#data = data.reshape(101, 61, 61)
data = data.reshape(61, 61, 101)
#data = data.reshape(61, 101, 61)
mat = data[:, 30, :]

fig, ax = plt.subplots()
ax.pcolor(mat, cmap = plt.cm.Blues)
plt.show()
