## Python helper file to plot and store graphs to represent scaling in distributed simulations
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

min_num_processes = 2
max_num_processes = 32
process_step = 3
show_error = False

folder = "sharc_timing_data/"
experiment = "w_exclfalse_60wait"

colours = ["blue","yellow","red"] # as many as the number of experiments included
legend_labels = ['Mean','Min','Max']
functions = [np.mean,np.min,np.max]

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [0,1,5,6]
scale_plots = [1,1,1,1,1,1]
plot_title = ["60 wait - excl=false. Simulation time","60 wait - excl=false. Computation time","60 wait - excl=false. Sync time","60 wait - excl=false. Idle time"]
plot_xlabel = ["Number of processes","Number of processes","Number of processes","Number of processes"]
plot_ylabel = ["Time (s)","Time (s)","Time (s)","Time (s)"]
plot_name = ["_sim_time.png","_comp_time.png","_sync_time.png","_idle_time.png"]

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

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y, error,title,xlabel,ylabel,name,colour,legend,show):
	if show_error:
		plt.errorbar(x, y, error,linewidth=2,color=colour,label=legend)
	else:
		plt.errorbar(x, y,linewidth=2,color=colour,label=legend)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.legend(loc=2)
	if show:
		plt.savefig(folder + name)
		plt.show()

# one plot per column per function
for i in range(len(columns_to_plot)):
	plt.figure()
	for f in range(len(functions)):
		# creating figure for the column plot
		
		w = (max_num_processes - min_num_processes+1) / process_step
		h = len(columns_to_plot)
		means = [[] for y in range(h)]
		stdevs = [[] for y in range(h)]
			
		for p in range(min_num_processes, max_num_processes+1,process_step):
			data = get_data_from_csv(folder + experiment + "__" + str(p))
			means[i].append(functions[f](data[:,columns_to_plot[i]]))
			stdevs[i].append(np.std(data[:,columns_to_plot[i]]))
		# means contains a list of arrays (one per columns analysed) with the average values per column across files
		# stdevs the same but with the corresponding st deviation

		#plot each column separately
		y = [x * scale_plots[i] for x in means[i]]
		error = [x * scale_plots[i] for x in stdevs[i]]
		plot(range(min_num_processes,max_num_processes+1,process_step),y,error,plot_title[i],plot_xlabel[i],plot_ylabel[i],experiment + plot_name[i],colours[f],legend_labels[f], f == len(functions)-1)	
		
		

    
