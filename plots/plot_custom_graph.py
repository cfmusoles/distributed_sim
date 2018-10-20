## Python helper file to plot and store graphs based on equations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

geometric_scaling = True
min_num_processes = 48
# for linear scaling of processors
max_num_processes = 288
process_step = 32
#for geometric scaling of processors
num_experiments = 7
geometric_step = 2

# simulation parameters
time_steps=7500 
total_neurons = 77000
average_neuron_connectivity = 3000
spike_frequencies = [20,80,160] # each one will generate one line

colours = ["red","green","blue","purple"] # as many as the number of spike_frequencies + 1
legend_labels = ['Handshake','Send / receive data (Low activity)','Send / receive data (Medium activity)','Send / receive data (High activity)']
plot_line_styles=['-','--','--','--','--']

# Each element on the following arrays corresponds to a column in columns_to_plot
plot_title = "Data exchange size per phase"
plot_xlabel = "Number of processes"
plot_ylabel = "Data sent"
plot_name = "a1"

# general plot settings
image_format = 'pdf'
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

def plot(x,y,linestyle,title,xlabel,ylabel,name,colour,legend,show):
	plt.errorbar(x, y,linestyle=linestyle,linewidth=2,color=colour,label=legend)
	if geometric_scaling:
		plt.yscale("linear",basey=10)
		plt.xscale("linear",basex=10)
	else:
		plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.legend(loc=2)
	if show:
		plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
		plt.show()

# creating figure for the column plot
plt.figure()
	
# X axis --> number of processes
if geometric_scaling:
	w = num_experiments
	x_axis = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	w = (max_num_processes - min_num_processes+1) / process_step
	x_axis = range(min_num_processes, max_num_processes+1,process_step)
	
# y axis --> cost of communication (handshake + send/receive)
y_axis = []
# handshake cost
y_axis.append([p*(p-1)*time_steps for p in x_axis ])
for freq in spike_frequencies:
	y_axis.append([p*(total_neurons / p * (freq * time_steps / 10000) * (p-1)) if total_neurons / p * average_neuron_connectivity > p else p*(total_neurons / p * (freq * time_steps / 10000) * total_neurons / p * average_neuron_connectivity) for p in x_axis ])


for i in range(len(y_axis)):
	plot(x_axis,y_axis[i],plot_line_styles[i],plot_title,plot_xlabel,plot_ylabel,plot_name,colours[i],legend_labels[i],i==len(y_axis)-1)
	


    
