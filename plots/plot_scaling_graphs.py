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

show_error = True
as_bar_plot = False

folder = "../results/archer/data/"
# each element on the following arrays corresponds to an experiment run (collection of files)
experiments = ["sparse_comm_st_randomBalanced_pex_pruned","sparse_comm_st_randomBalanced_nbx_pruned","sparse_comm_st_hypergraphPartitioning_pex_pruned","sparse_comm_st_hypergraphPartitioning_nbx_pruned"] # plot more than one set of results in the graphs
colours = ["red","green","blue","orange"] # as many as the number of experiments included
legend_labels = ['Random-PEX','Random-NBX','Hypergraph partition-PEX','Hypergraph partition-NBX']

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [1,2,3,6,7]
reference_values = [1,2,1,6,7] # used to take values on each column divided by these
use_ref_values = False
scale_plots = [1,1,1,1,1,1]
plot_title = ["Simulation time","Computation time","Computation variance","Data Exchange time","Implicit sync time"]
plot_xlabel = ["Number of processes","Number of processes","Number of processes","Number of processes","Number of processes","Number of processes"]
plot_ylabel = ["Time(s)","Time(s)","Variance","Time(s)","Time(s)","Time(s)"]
image_format = 'pdf'
plot_name = ["a1","a2","a3","a4","a5"]

bar_plot_size = 0.5 / len(experiments)

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

def plot(x,y, error,title,xlabel,ylabel,name,colour,legend,show,global_counter):
	if as_bar_plot:
		rx = x + global_counter * bar_plot_size*np.array(x)
		plt.bar(rx,y,width=bar_plot_size*np.array(x),color=colour,label=legend)
		
	else:
		if show_error:
			plt.errorbar(x, y, error,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
		else:
			plt.errorbar(x, y,linewidth=1,color=colour,label=legend,marker='s',markersize=5)
	if geometric_scaling:
		plt.yscale("log",basey=10)
		plt.xscale("log",basex=10)
	else:
		plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.xticks(experiment_range,experiment_range)
	plt.tick_params(axis='x',which='minor',bottom=False)
	#plt.tight_layout()
	plt.gcf().subplots_adjust(left=0.17)
	if len(experiments) > 1:
		plt.legend(loc='best')
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
			data = get_data_from_csv(folder + experiments[j] + "__" + str(p))
			if use_ref_values:
				ref = [row[columns_to_plot[i]] / row[reference_values[i]] for row in data]
				means[i].append(np.mean(ref))
				stdevs[i].append(np.std(ref))
			else:
				means[i].append(np.mean(data[:,columns_to_plot[i]]))
				stdevs[i].append(np.std(data[:,columns_to_plot[i]]))
		# means contains a list of arrays (one per columns analysed) with the average values per column across files
		# stdevs the same but with the corresponding st deviation

		#plot each column separately
		y = [x * scale_plots[i] for x in means[i]]
		error = [x * scale_plots[i] for x in stdevs[i]]
		plot(experiment_range,y,error,plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colours[j],legend_labels[j], j == (len(experiments)-1),j)	


    
