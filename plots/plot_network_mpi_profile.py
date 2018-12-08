## Python helper file to plot and store graphs to represent network bandwidth and MPI performance

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from pylab import *

num_processes = 192
plot_bandwidth = True			# plot network bandwidth data
plot_sent_data = False			# plot application sent data
plot_comm_cost = False			# plot combined comm cost
storeResults = False

folder = "../../hyperPraw/"
bandwidth_send_experiment_name = 'results_mpi_send_bandwidth_' + str(num_processes)
sim_sent_experiment = 'without_prediction_comm_matrix_' + str(num_processes)

xlabel = "Process"
ylabel = "Process"
zlabels = ["MB/sec","Bytes"]
titles = ["P2P Bandwidth (praw-noB)","Sent data (praw-withB)","Cost of communication (praw-withB)"]
filenames = ["a1","a2","a3"]
image_format = 'pdf'

# general plot settings
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'small',
        'legend.fontsize': 'small',
        'font.size': 14,
        'savefig.dpi': 200,
}
plt.rcParams.update(fig_settings)

def get_data_from_csv(filename,delimiter=","):
	data = np.genfromtxt(filename,skip_header=0,delimiter=delimiter)
	return data

def plot_3dgraph(data,title,zlabel,filename):
	fig = plt.figure()
	ax = Axes3D(fig)

	lx= len(data[0])            # Work out matrix dimensions
	ly= len(data[:,0])
	xpos = np.arange(0,lx,1)    # Set up a mesh of positions
	ypos = np.arange(0,ly,1)
	xpos, ypos = np.meshgrid(xpos+0.5, ypos+0.5)

	xpos = xpos.flatten()   # Convert positions to 1D array
	ypos = ypos.flatten()
	zpos = np.ones(lx*ly)*1e-10

	dx = 1. * np.ones_like(zpos)
	dy = dx.copy()
	dz = data.flatten()

	colors = plt.cm.jet(data.flatten()/float(data.max()))
	ax.bar3d(xpos,ypos,zpos, dx, dy, dz,  color=colors, alpha=1., zsort='max',linewidth=0)
	
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_zlabel(zlabel)
	
	plt.title(title)
	if storeResults:
		plt.savefig(filename + "." + image_format,format=image_format,dpi=300)
	plt.show()



## creating figure for bandwidth estimation (send performance)
if plot_bandwidth:
	data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	plot_3dgraph(data,titles[0],zlabels[0],filenames[0])

## creating combined figure for cost of communication
## each value corresponds to a process-process pair: total data sent / bandwidth 
if plot_sent_data:
	data = get_data_from_csv(folder + sim_sent_experiment,' ')
	plot_3dgraph(data,titles[1],zlabels[1],filenames[1])

if plot_comm_cost:
	bandwidth_data = get_data_from_csv(folder + bandwidth_send_experiment_name,'\t')
	comm_data = get_data_from_csv(folder + sim_sent_experiment,' ')
	cost_data = np.array([c/b for b, c in zip(bandwidth_data, comm_data)])
	cost_data = np.nan_to_num(cost_data)
	print('Total cost: ' + str(cost_data.sum()))
	# colour maps http://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html
	pcolor(cost_data,cmap=get_cmap("binary"))
	cbar = colorbar()
	plt.title(titles[2])
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	cbar.ax.set_ylabel("Cost (time)",rotation=90)
	if storeResults:
		plt.savefig(filenames[2] + "." + image_format,format=image_format,dpi=300)
	plt.show()

