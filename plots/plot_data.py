## Python helper file to plot and store graphs to represent neuronal activity
# Graphs produced:
# - Total & per population number of spikes histogram
# - Total & per population Inter-spike interval histogram
# - ISI coefficient of variation across populations (http://neuronaldynamics.epfl.ch/online/Ch7.S3.html)
# - Synchrony across populations (Potjans and Diesmann, 2014)

# TODO
# -scale in Brette paper for ISI is 3...10...30...100...300...1000 (x[3 + 1/3])

import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec

experiment = "cm"
populations =  ["l2_3e","l2_3i","l4e","l4i","l5e","l5i","l6e","l6i"] # ["Exc","Inh"]
num_processes = 4
sim_time = 0.5
force_plot_format = False	# plot formatting matches VAbenchmarks if true

ISI_filename = experiment + "__" + str(num_processes) + "_ISI"
spikes_filename = experiment + "__" + str(num_processes) + "_spikes"
isicv_filename = experiment + "__" + str(num_processes) + "_ISICV"

# general plot settings
image_format = 'pdf'
plt.rcParams['figure.facecolor'] = 'white'

def plot_ISI_histogram(data,name):
	if(len(data) == 0):
		 return
	if force_plot_format:
		#bins = np.logspace(math.log10(1), math.log10(3000),base=10) 
		bins_log = np.arange(0, 8, 0.2)
		bins = np.exp(bins_log)
		if name == "Exc":
			#plt.xticks(np.log([10, 100, 1000]))
			plt.yticks(np.arange(0,2000,200))
			plt.ylim(0, 1800)
			plt.xlim(1,3000)
			plt.text(900, 1600, 'exc')
		if name == "Inh":
			#plt.xticks=np.log([10, 100, 1000])
			plt.yticks(np.arange(0,500,50))	
			plt.ylim(0, 450)
			plt.xlim(1,1500)
			plt.text(900, 350, 'inh')		
	else:
		bins = np.logspace(math.log10(min(data)), math.log10(max(data)),base=10)
 
	plt.hist(data,bins=bins) 
	plt.xscale('log')
	plt.ylabel('n in bin')
	if not force_plot_format:
		plt.title(name + ': Inter-spike Interval graph')
		plt.xlabel('Inter-spike interval (ms)')
	plt.grid(False)
	# plt.axis([40, 160, 0, 0.03]) # to define x and y min max values
	ax = plt.axes()
	ax.xaxis.set_major_formatter(ScalarFormatter())
	plt.savefig("ISI_" + name + "." + image_format,format=image_format,dpi=1000)
	plt.show()

def plot_spikes_histogram(data,name):
	if(len(data) == 0):
		return
	binwidth = 1
	plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
	plt.xlabel('Number of spikes')
	plt.ylabel('N in bin')
	plt.title(name + ': Spikes graph')
	plt.grid(False)
	# plt.axis([40, 160, 0, 0.03]) # to define x and y min max values
	plt.savefig("spikes_" + name + "." + image_format,format=image_format,dpi=1000)
	plt.show()

def plot_ISICV_histogram(data,name):
	if(len(data) == 0):
		return
	if force_plot_format:
		bin_width = 0.1
		bins = np.arange(0, 2, bin_width)
		xticks = np.arange(0, 2, 0.5)
		plt.xticks(xticks)
		plt.hist(data, bins=bins)
		if name == "Exc":
			plt.yticks(np.arange(0,500,50))
			plt.ylim(0, 450)
			plt.text(1.1, 350, 'exc')
		if name == "Inh":
			plt.yticks(np.arange(0,140,20))
			plt.ylim(0, 120)
			plt.text(1.1, 100, 'inh')
	else:
		bin_width = 0.1
		plt.hist(data, bins=np.arange(min(data), max(data) + bin_width, bin_width))
	plt.ylabel('n in bin')
	if not force_plot_format:
		plt.title(name + ': ISI CV graph')
		plt.xlabel('CV(ISI)')
	plt.grid(False)
	# plt.axis([40, 160, 0, 0.03]) # to define x and y min max values
	plt.savefig("ISICV_" + name + "." + image_format,format=image_format,dpi=1000)
	plt.show()

### all neurons plot
plot_ISI_histogram(np.loadtxt(ISI_filename),"Total")
plot_spikes_histogram(np.loadtxt(spikes_filename),"Total")
plot_ISICV_histogram(np.loadtxt(isicv_filename),"Total")

### ISI, spikes and ISI-CV per population
spike_pops = []
isi_pops = []
isicv_pops = []
for pop in populations:
	spikedata = np.loadtxt(spikes_filename + "_" + pop)
	spike_pops.append(spikedata) #[x/sim_time for x in spikedata])
	isidata = np.loadtxt(ISI_filename + "_" + pop)
	isi_pops.append(isidata)
	isicvdata = np.loadtxt(isicv_filename + "_" + pop)
	isicv_pops.append(isicvdata)
	plot_spikes_histogram(spikedata,pop)
	plot_ISI_histogram(isidata,pop)
	plot_ISICV_histogram(isicvdata,pop)

### all populations freq spikes

fig = plt.figure()
ax = plt.axes()
#plt.hold(True)
counter = 1
for pop in spike_pops:
	if(len(pop) > 0):
		pop = [x/sim_time for x in pop]
		bp = plt.boxplot(pop, positions=[counter], widths = 0.6)
		elements = ['boxes','caps','whiskers']
		color = 'blue'
		# Iterate over each of the elements changing the color
		for elem in elements:
			[plt.setp(bp[elem][idx], color=color) for idx in xrange(len(bp[elem]))]
		counter += 1

ax.set_xticklabels([""] + populations + [""])
ax.set_xticks(range(0,len(populations)+2))
plt.xlabel('Population')
plt.ylabel('Freq of spikes (Hz)')
plt.title('Frequency of spikes per population')
plt.savefig(experiment + "_freq." + image_format,format=image_format,dpi=1000)
plt.show()

### all populations synchrony (on spike count histogram, variance over mean)
y_pos = np.arange(len(populations))
synchrony = []
for pop in spike_pops:
	synchrony.append(np.var(pop) / np.mean(pop))
bar_width = 0.4 
plt.bar(y_pos, synchrony, bar_width, align='center', alpha=0.5)
plt.xticks(y_pos, populations)
plt.ylabel('Synchrony (variance / mean)')
plt.title('Populations')
plt.savefig(experiment + "_synchrony." + image_format,format=image_format,dpi=1000)
plt.show()

### all populations coefficient of variation (on the ISI, std over the mean)
y_pos = np.arange(len(populations))
cv = []
for pop in isi_pops:
	cv.append(np.std(pop) / np.mean(pop))
bar_width = 0.4 
plt.bar(y_pos, cv, bar_width, align='center', alpha=0.5)
plt.xticks(y_pos, populations)
plt.ylabel('ISI CV (std / mean)')
plt.title('Populations')
plt.savefig(experiment + "_ISICV." + image_format,format=image_format,dpi=1000)
plt.show()




