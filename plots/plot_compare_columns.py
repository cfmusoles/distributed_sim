## Python helper file to plot and store graphs comparing the results from two experiments on two columns
# Graphs produced: produces a two column comparison on the differences on those two columns for the two experiments

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import seaborn as sns

geometric_scaling = True
min_num_processes = 192
# for linear scaling of processors
max_num_processes = 288
process_step = 32
#for geometric scaling of processors
num_experiments = 6
geometric_step = 2

show_title = False

folder = "../results/archer/mvc160/"
# each element on the following arrays corresponds to an experiment run (collection of files)
# only 2 are acceptable
experiments = ["mvc160_roundrobin_pex","mvc160_hypergraphPartitioning_nbx"]

# Each element on the following arrays corresponds to a column in columns_to_plot
columns_to_plot = [0,1]
columns_names = ['Build cost','Sim time gain']
colours = ["red","green"]
scale_plots = [0,0,0,1,0] # if 0, show difference in percentage
plot_title = "Build time cost (HB-NBX over round robin)"
plot_xlabel = "Number of processes"
plot_ylabel = "Time (s)"
image_format = 'pdf'
plot_name = 'a1'

if geometric_scaling:
	experiment_range = [min_num_processes * geometric_step ** (n-1) for n in range (1, num_experiments+1)]
else:
	experiment_range = range(min_num_processes, max_num_processes+1,process_step)

# load data from both columns
def get_data_from_csv(filename):
	data = np.genfromtxt(filename,skip_header=1,delimiter=",")
	return data

h = len(columns_to_plot)
means = [[] for y in range(h)]

# one plot per column
for i in range(len(columns_to_plot)):	
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
		comp = data1[:,columns_to_plot[i]] - data2[:,columns_to_plot[i]]
		means[i].append(np.mean( np.abs(comp) ))


# construct data structure and plot
df = pd.DataFrame({
    'Series': experiment_range,
    columns_names[0]: means[0],
    columns_names[1]: means[1]
})
fig, ax1 = plt.subplots()
tidy = df.melt(id_vars='Series').rename(columns=str.title)
sns.barplot(x='Series', y='Value', hue='Variable', data=tidy, ax=ax1,palette=colours)
ax1.set(xlabel=plot_xlabel, ylabel=plot_ylabel)
l = ax1.legend()
l.set_title('')
plt.savefig(plot_name + "." + image_format,format=image_format,dpi=1000)
plt.show()