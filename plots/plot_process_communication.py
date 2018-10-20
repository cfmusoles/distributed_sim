# Plots the Process Communication graph for a simulation result
# Shows the frequency of communication of all edges

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import networkx as nx

num_processes = 8

folder = "../"
experiment_name = "test" 
row_selected = 0
column_selected = 6
plot_name = "partition_communication_nonempty_runtime.png"

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

def get_data_from_csv(filename,skip_header=1,delimiter=","):
	data = np.genfromtxt(filename,skip_header=skip_header,delimiter=delimiter)
	return data

def get_data_from_csv_irregular_columns(filename,delimiter=",",skip_header=True):
    import csv
    datafile = open(filename, 'r')
    datareader = csv.reader(datafile)
    if(skip_header):
        next(datareader,None)
    data = []
    for row in datareader:
        d = [ elem for elem in row[0].split(delimiter) ]
        if '' in d:
            d.remove('')
        d = [int(elem) for elem in d]
        data.append( d )
    return data

# load hyperedge strengths (based on runtime frequency of communication)
hyperedges_strength = []
for p in range(num_processes):
    filename = folder + experiment_name + "_" + str(num_processes) + "_" + str(p)
    data = get_data_from_csv(filename)
    if len(data.shape) > 1:
        hyperedges_strength.append(data[row_selected][column_selected])
    else:
        hyperedges_strength.append(data[column_selected])

# load processor graph (based on process connectivity) and create NX graph
G = nx.DiGraph()
filename = folder + experiment_name + "_partition_graph_" + str(num_processes)
connectivity = get_data_from_csv_irregular_columns(filename,delimiter=" ")
hyperedge = []
for n in range(len(connectivity)):
    G.add_node(n)
    hyperedge = []
    for c in range(1,len(connectivity[n])):
        hyperedge.append((n,connectivity[n][c]))
    G.add_edges_from(hyperedge,color=hyperedges_strength[n])

edge_colours = [(0.5,1-G[u][v]['color'],0,1) for u,v in G.edges()]
nx.draw_networkx(G,with_labels=True,node_color='red',node_size=500,edge_color=edge_colours,edge_cmap='viridis')


#color bar legend
#fig = plt.gcf()
#ax = plt.gca()
fig, ax = plt.subplots()
cmap = mpl.colors.LinearSegmentedColormap.from_list("MyCmapName",[(0.5,1,0,1),(0.5,0,0,1)])
norm = mpl.colors.Normalize(vmin=0, vmax=100)
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Frequency of communication(%)')
fig.show()

plt.draw()
plt.show()
#plt.savefig(folder + name)