# Plots bar graph with proportional results (relative to a value) in scaling results (increasing number of processes)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

geometric_scaling = True
min_num_processes = 192
# for linear scaling of processors
max_num_processes = 32
process_step = 3
#for geometric scaling of processors
num_experiments = 6
geometric_step = 2

folder = "../results/archer/mvc160/"
experiment = "mvc160_hypergraphPartitioning_nbx" #"mvc160_hypergraphPartitioning_nbx"

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [1,2,6,7] #reference column is the first one
useReferenceValue = False
expressAsPercentage = True
show_plot_title = False
colours = 'rbc'#'byrwgcm' 
legend_labels = ['Computation','Data Exchange','Implicit sync']

plot_title = "Sim time (HB-NBX)"
plot_xlabel = "Number of processes"
plot_ylabel = "Time (% of total)"
plot_name = "a1"

# general plot settings
image_format = 'pdf'
plt.rcParams['figure.facecolor'] = 'white'
fig_settings = {  
        'lines.linewidth': 0.5,
        'axes.linewidth': 0.5,
        'axes.labelsize': 'large',
        'legend.fontsize': 'medium',
        'font.size': 16,
        'savefig.dpi': 1000,
}
plt.rcParams.update(fig_settings)

def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

def plot(x,y):
	x_pos = np.arange(len(x))

	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	patch_handles = []
	bottom = np.zeros(len(x)) # left alignment of data starts at zero
	for i, d in enumerate(y):
		patch_handles.append(ax.bar(x_pos, d, 
			color=colours[i%len(colours)], align='center', 
			bottom=bottom,label=legend_labels[i%len(legend_labels)]))
		# accumulate the left-hand offsets
		bottom += d
	
	ax.set_xticks(x_pos)
	ax.set_xticklabels(x)
	plt.xlabel(plot_xlabel)
	plt.ylabel(plot_ylabel)
	if show_plot_title:
		plt.title(plot_title)
	#plt.legend(loc=3)
	plt.legend(bbox_to_anchor=(1.0,1.15))
	plt.savefig(plot_name + "." + image_format,format=image_format,dpi=1000)
	plt.show()


w = (max_num_processes - min_num_processes+1) / process_step
h = len(columns_to_plot)
means = [[] for y in range(h)]

if geometric_scaling:
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)
	
for i in range(len(columns_to_plot)):
	# creating figure for the column plot
	for p in experiment_range:
		data = get_data_from_csv(folder + experiment + "__" + str(p))
		means[i].append(np.mean(data[:,columns_to_plot[i]]))
		# means contains a list of arrays (one per columns analysed) with the average values per column across files
		
# convert column average values to proportions relative to first column reference
if useReferenceValue:
	for i in range(1,len(means)):
		for w in range(len(means[i])):
			means[i][w] = means[i][w] / means[0][w] * 100
a = means.pop(0) # no longer need reference values
if expressAsPercentage:
	for row,mean in enumerate(means):
		means[row] = [e / a[i] * 100 for i,e in enumerate(mean)]


#plot bar graph (all processor counts and columns)
plot(experiment_range,means)	


    
