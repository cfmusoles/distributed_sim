import numpy as np
import matplotlib.pyplot as plt
import math

###### SETUP######
experiment_folder = "../"
experiment_names = ["test"]
total_processes = 12
row_selected = 1	#if you have more than one row per processor (-1 for average)
isTrace = True		# if true, plot will represent real processor traces; false for average plots
# colours correspond to each column being considered
if isTrace:
	colors ='rrrrybr'#'rybwgcm' #different functions being tracked
else:
	colors ='rybwgcm' # elements timed (comp, idle, sync...)
## for tracing timings ##
#how many measurements (from start to end cycles (for loops in sim)) 
loop_size = 7 # how many measurements are traced (from distSim.cpp)
cycle_size = 1 # to see sync events, make this the propagation_time_step (distSim.cpp)
start_cycle = 0 #inclusive (multiple of cycle_size to see send to send events)
num_cycles = 10
## for average timings ##
#how many measurements (compute, idle and sync); column selection (inclusive)
start_col = 0
end_col = 2
#################

image_format = 'pdf'


if isTrace:
	start_column = start_cycle * (loop_size * cycle_size)
	end_column = start_column + (loop_size * cycle_size * num_cycles) - 1
else:
	start_column = start_col
	end_column = end_col

def get_data_from_csv(filename):
	range_value = 0
	data = [[] for _ in range(start_column,end_column+1)]
	for p in range(0,total_processes):
		if isTrace:
			mydata = np.genfromtxt(filename + "_" + str(p) + "_trace",skip_header=0,delimiter=",")
		else:
			mydata = np.genfromtxt(filename + "_" + str(p),skip_header=1,delimiter=",")
		mydata = np.array(mydata)
		if len(mydata.shape) > 1:
			mydata = [[each_list[i] for i in range(start_column,end_column+1)] for each_list in mydata]
		else:
			mydata = [mydata[i] for i in range(start_column,end_column+1)]
		mydata = np.array(mydata)
		if row_selected == -1:
			mydata = np.mean(mydata,axis=0)
		else:
			if len(mydata.shape) > 1:
				mydata = mydata[row_selected]
		if isTrace: #convert data from timeline to time spent
			mydata = [(mydata[i]-mydata[i-1]) for i in range(1,len(mydata))]
			mydata.append(0)
		m = sum(mydata)
		if m > range_value:
			range_value = m
		for c in range(0,len(mydata)):
			data[c].append(mydata[c])
	
	return data, range_value

datas = []
max_value = 0
for e in range(0,len(experiment_names)):
	data, m = get_data_from_csv(experiment_folder + experiment_names[e] + "_" + str(total_processes))
	datas.append(data)
	if m > max_value:
		max_value = m

for ds in range(0,len(datas)):
	processes = range(0,total_processes)
	# generate some multi-dimensional data & arbitrary labels
	#percentages = (np.random.randint(5,20, (len(processes), number_of_timings)))
	y_pos = np.arange(len(processes))

	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111)

	patch_handles = []
	left = np.zeros(len(processes)) # left alignment of data starts at zero
	for i, d in enumerate(datas[ds]):
		patch_handles.append(ax.barh(y_pos, d, 
			color=colors[i%len(colors)], align='center', 
			left=left))
		# accumulate the left-hand offsets
		left += d

	# go through all of the bar segments and annotate
	#for j in xrange(len(patch_handles)):
	#	for i, patch in enumerate(patch_handles[j].get_children()):
	#		bl = patch.get_xy()
	#		x = 0.5*patch.get_width() + bl[0]
	#		y = 0.5*patch.get_height() + bl[1]
	#		ax.text(x,y, "%ds" % (percentages[i][j]), ha='center')

	ax.set_yticks(y_pos)
	ax.set_yticklabels(processes)
	ax.set_xlabel('Time (s)')
	ax.set_ylabel("Process ID")
	#if max_value < 1:
	#	top_value = max_value * 1.5
	#	xlabel_tick_step = max_value/len(datas[ds])
	#else:
	#	xlabel_tick_step = math.floor(max_value/len(datas[ds]) * 10)/10
	#	top_value = math.ceil(max_value) + xlabel_tick_step
	#	if not top_value % xlabel_tick_step == 0:
	#		top_value += top_value % xlabel_tick_step
		
	#ax.set_xticks(np.arange(0,top_value,xlabel_tick_step))
	
	plt.title(experiment_names[ds])
	plt.savefig(experiment_names[ds] + "." + image_format,format=image_format,dpi=1000)
	plt.show()
