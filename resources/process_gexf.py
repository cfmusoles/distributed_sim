#install networkx: pip install --user networkx

import networkx as nx
import matplotlib.pyplot as plt
import math
import os
import random

graph_path = "v1_2_reduced.gexf"
processed_path = "v1_2_reduced.txt"
inh_exc_ratio = 0.1

graph = nx.read_gexf(graph_path)

print nx.info(graph)

graph = nx.convert_node_labels_to_integers(graph,0)
nx.write_adjlist(graph,processed_path)

def prepend_and_filter_file(originalfile,string):
	f2 = open('newfile.txt','w')
	with open(originalfile,'r') as f:
		f2.write(string)
		for line in f:
			if '#' not in line:
				f2.write(line)
	os.rename('newfile.txt',originalfile)


# remove comment lines at the top and add header information
number_nodes = nx.number_of_nodes(graph)
num_inh = int(math.floor(number_nodes * inh_exc_ratio))
num_exc = number_nodes - num_inh
prepend_and_filter_file(processed_path, str(number_nodes) + " " + str(num_exc) + " " + str(-num_inh) + "\n")

# assign neurons to pop1 (exc) or pop2 (inh) randomly
inh_cells = [str(random.randint(0,number_nodes-1)) for _ in range(num_inh)]
temp = open('temp.txt','w')
f = open(processed_path,'r')
counter = 0
for line in f:
	finalLine = ""
	if counter == 0:
		# leave header as is
		finalLine = line
	else:
		tokens = line.split()
		if tokens[0] in inh_cells:
			finalLine += "1"
		else:
			finalLine += "0"
		for t in range(1,len(tokens)):
			finalLine += " " + tokens[t]
		finalLine += "\n"
	temp.write(finalLine)
	counter+=1
		
temp.close()
f.close()
os.rename('temp.txt',processed_path)

