## Python helper file to plot and store graphs to represent scaling in distributed simulations
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

geometric_scaling = True
min_num_processes = 96
# for linear scaling of processors
max_num_processes = 288
process_step = 32
#for geometric scaling of processors
num_experiments = 6
geometric_step = 2

row_selected = 0 #for data files with more than one experiment, select one row

show_error = True

folder = "../results/archer/data/"
# each element on the following arrays corresponds to an experiment run (collection of files)
experiments = ["sparse_comm_st_randomBalanced_pex","sparse_comm_st_randomBalanced_nbx","sparse_comm_st_hypergraphPartitioning_pex","sparse_comm_st_hypergraphPartitioning_nbx"] # plot more than one set of results in the graphs
colours = ["red","green","blue","yellow"] # as many as the number of experiments included
legend_labels = ['Random Balanced-PEX','Random Balanced-NBX','Hypergraph partition-PEX','Hypergraph partition-NBX']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [0,1,2]
reference_columns = [0,2]
scale_plots = [1,1,1]
plot_title = ["Computational time","Implicit sync time","Data exchange time"]
plot_xlabel = ["Number of processes","Number of processes","Number of processes"]
plot_ylabel = ["Time(s)","Time(s)","Time(s)"]
image_format = 'pdf'
plot_name = ["a1","a2","a3"]

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

if geometric_scaling:
	w = num_experiments
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	w = (max_num_processes - min_num_processes+1) / process_step
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, error,title,xlabel,ylabel,name,colour,legend,show):
	if show_error:
		#plt.errorbar(x, y, error,linewidth=0,color=colour,label=legend,marker='s',markersize=5)
		plt.errorbar(x,y,yerr=error,fmt='s-',label=legend,color=colour)
	else:
		plt.errorbar(x,y,fmt='s',label=legend,color=colour)
	if geometric_scaling:
		plt.yscale("log",basey=10)
		plt.xscale("log",basex=10)
	else:
		plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks(experiment_range,experiment_range)
	plt.tick_params(axis='x',which='minor',bottom=False)
	plt.title(title)
	if len(experiments) > 1:
		plt.legend(loc=2)
	if show:
		plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
		plt.show()

# one plot per column
for i in range(len(columns_to_plot)):
	# creating figure for the column plot
	plt.figure()
	for j in range(len(experiments)):
		h = len(columns_to_plot)
		means = [[] for y in range(h)]
		stdevs = [[] for y in range(h)]
		
		for p in experiment_range:
			experiment_data = []
			for n in range(p):
				data = get_data_from_csv(folder + experiments[j] + "_" + str(p) + "_" + str(n))
				experiment_data.append(data[row_selected,columns_to_plot[i]])
				
			means[i].append(np.mean(experiment_data))
			stdevs[i].append(np.std(experiment_data))
		# means contains a list of arrays (one per columns analysed) with the average values per column across files
		# stdevs the same but with the corresponding st deviation

		#plot each column separately
		y = [x * scale_plots[i] for x in means[i]]
		error = [x * scale_plots[i] for x in stdevs[i]]
		plot(experiment_range,y,error,plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colours[j],legend_labels[j], j == (len(experiments)-1))	


    
