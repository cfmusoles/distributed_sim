## Python helper file to plot and store graphs comparing the results from two experiments
# Graphs produced: as many as columns that want to be compared

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

geometric_scaling = True
min_num_processes = 192
# for linear scaling of processors
max_num_processes = 288
process_step = 32
#for geometric scaling of processors
num_experiments = 6
geometric_step = 2
colour = "green"

show_title = False

folder = "../results/archer/mvc160/"
# each element on the following arrays corresponds to an experiment run (collection of files)
# only 2 are acceptable
experiments = ["mvc160_roundrobin_pex","mvc160_hypergraphPartitioning_nbx"] # plot more than one set of results in the graphs
legend_label = ""#HP gain over Random"

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [16,8,22,0,1] # 8, 3, 22
scale_plots = [0,0,0,1,0] # if 0, show difference in percentage
reference_values = [17,8,23,0,0] # used to take values on each column divided by these
use_ref_values = True
plot_title = ["Remote spikes (HB-NBX over round robin)","Data volume (HP-NBX over Round robin)","ARN reduction (Round Robin vs hypergraph)","Build time cost (HB-NBX over round robin)","Simulation time gain (HB-NBX over round robin)"]
plot_xlabel = ["Number of processes","Number of processes","Number of processes","Number of processes","Number of processes"]
plot_ylabel = ["Spikes sent difference (%)","Data exchanged difference (%)","ARN difference (%)","Time loss (s)","Time gain (s)"]
image_format = 'pdf'
plot_name = ['a' + str(x) for x in range(len(columns_to_plot))] #["a1","a2","a3","a4","a5","a6","a7"]


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

def plot(x,y, error,title,xlabel,ylabel,name,colour,legend):

	if geometric_scaling:
		plt.bar(x,y,width=np.array(x)/2,color=colour,label=legend,align="edge")
		plt.yscale("linear",basey=10)
		plt.xscale("log",basex=10)
	else:
		plt.bar(x,y,width=10.55,color=colour,label=legend)
		plt.yscale("linear")
		plt.xscale("linear")
	#plt.errorbar(x, y,linewidth=2,color=colour,label=legend)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.xticks([e*1.25 for e in experiment_range],experiment_range)
	plt.tick_params(axis='x',which='minor',bottom=False,labelbottom=False)
	if show_title:
		plt.title(title)
	#plt.ylim(0,100)
	plt.gcf().subplots_adjust(left=0.18,bottom=0.2)
	if not legend == "":
		plt.legend(loc=1)
	plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
	plt.show()


if geometric_scaling:
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)

# one plot per column
for i in range(len(columns_to_plot)):
	# creating figure for the column plot
	plt.figure()
	h = len(columns_to_plot)
	means = [[] for y in range(h)]
	
	for p in experiment_range:
		data1 = get_data_from_csv(folder + experiments[0] + "__" + str(p))
		data2 = get_data_from_csv(folder + experiments[1] + "__" + str(p))
		#make sure both sets of experiments contain the same number of repetitions (drop odd ones)
		rem = abs(len(data1) - len(data2))
		if rem > 0:
			direction = len(data1) > len(data2)
			for r in range(rem):
				if direction:
					data1 = np.delete(data1,(len(data1)-1),axis=0)
				else:
					data2 = np.delete(data2,(len(data2)-1),axis=0)
		if scale_plots[i] == 0:
			comp = (data1[:,columns_to_plot[i]] - data2[:,columns_to_plot[i]]) / data1[:,columns_to_plot[i]] * 100
		else:
			comp = data1[:,columns_to_plot[i]] - data2[:,columns_to_plot[i]]
		if use_ref_values:
			comp = comp / ( data2[:,reference_values[i]] - data1[:,reference_values[i]] ) * 100
		means[i].append(np.mean( comp ))
	
	#plot each column separately
	if scale_plots[i] == 0:
		y = means[i]
	else:	
		y = [x * scale_plots[i] for x in means[i]]
	plot(experiment_range,y,[],plot_title[i],plot_xlabel[i],plot_ylabel[i],plot_name[i],colour,legend_label)	
	

    
