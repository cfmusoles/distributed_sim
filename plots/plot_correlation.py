## Python helper file to plot and store graphs correlating two variables
## The first variable can be baselined (measuring improvement)

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

folder = "../results/archer/mvc160/"
# first experiment is the baseline
experiments = ["mvc160_roundrobin_pex","mvc160_hypergraphPartitioning_nbx"]
legend_label = "NBX gain over PEX"

# Use the baseline column (y axis) to plot improvements on compared column
baseline_column = 0
compared_column = 1
plot_title = "Correlation between runtime neighbours and sync time"
plot_xlabel = "Difference in runtime neighbours"
plot_ylabel = "Sync time improvement"
image_format = 'pdf'
plot_name = "a1"

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
	#plt.errorbar(x, y,linewidth=2,color=colour,label=legend)
	plt.scatter(x,y,marker='^',color=colour,label=legend)
	if geometric_scaling:
		plt.yscale("linear",basey=10)
		plt.xscale("linear",basex=10)
	else:
		plt.yscale("linear")
		plt.xscale("linear")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	if len(experiments) > 1:
		plt.legend(loc=1)
	plt.savefig(name + "." + image_format,format=image_format,dpi=1000)
	plt.show()


if geometric_scaling:
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)

# creating figure for the column plot
plt.figure()
baseline_values = []
compared_values = []

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
	baseline_values = (data1[:,baseline_column] - data2[:,baseline_column]) / data1[:,baseline_column] * 100
	compared_values = (data1[:,compared_column] - data2[:,compared_column]) / data1[:,compared_column] * 100

plot(baseline_values,compared_values,[],plot_title,plot_xlabel,plot_ylabel,plot_name,colour,legend_label)	


    
