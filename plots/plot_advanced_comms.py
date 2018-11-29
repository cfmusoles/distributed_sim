import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec

folder = "../"
experiments = ["test_comm_sizes"]
colours= ['blue','yellow']
legends = ['Random','Partition']
num_processes = 4

# general plot settings
image_format = 'pdf'
plt.rcParams['figure.facecolor'] = 'white'

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=0,delimiter=" ")
	return data


def plot_histogram(datas):
	for d in range(len(datas)):
		bins = range(int(math.floor(min(datas[d]))), int(math.ceil(max(datas[d]))),1)
	 	plt.hist(datas[d],bins=bins,color=colours[d],alpha=0.5,label=legends[d]) 
	#plt.xscale('log')
	plt.ylabel('n in bin')
	plt.title('Message sizes during simulation')
	plt.xlabel('Message size (B)')
	plt.grid(False)
	# plt.axis([40, 160, 0, 0.03]) # to define x and y min max values
	ax = plt.axes()
	ax.xaxis.set_major_formatter(ScalarFormatter())
	plt.legend(loc=1)
	plt.savefig("message_sizes" + str(num_processes) + "." + image_format,format=image_format,dpi=1000)
	plt.show()


### all message sizes plot
datas = []
for e in experiments:
	datas.append(get_data_from_csv(folder + e + "_" + str(num_processes)))

plot_histogram(datas)

