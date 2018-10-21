import csv
import networkx as nx
import os
from os import path
import numpy as np
import scipy.io
#import ipdb
from copy import deepcopy
from tqdm import tqdm, trange

    # Create version 1.1 graph file from the csv

def main():



    graph = 'v1_2.gexf'

    print "Loading Full Graph"
    G = nx.read_gexf(graph)
    #A = np.transpose(nx.adjacency_matrix(G,nodelist=nodes,weight=key).toarray())

    # get a list of all LPUs
    lpus = nx.read_gexf('v1_2_lpu.gexf').nodes();

    # for every LPU:
    #for lpu in lpus:

    for lpu in tqdm(lpus, total=len(lpus)):
        # create a seperate graph
        G_temp = deepcopy(G)
        #ipdb.set_trace()
        # create a new adjacency matrix
        rm_nodes = []
        for node in G_temp:
                
            if G_temp.node[node]['innv_neuropil'] != lpu:
                rm_nodes.append(node) 
        for node in rm_nodes:      
            G_temp.remove_node(node)

        Gc = max(nx.connected_component_subgraphs(G_temp.to_undirected()), key=len)
        removed = set(G_temp.nodes()) - set(Gc.nodes())
        for n in removed:
            G_temp.remove_node(n)

        A = (nx.adjacency_matrix(G_temp,nodelist=G_temp.nodes(),weight='weight').toarray())
        obj_arr = np.zeros((len(G_temp.nodes()),), dtype=np.object)
                       
        for i,n in enumerate(G_temp.nodes()):
            obj_arr[i] = n    
        scipy.io.savemat('lpus/v1_2_lpu_%s.mat' % lpu, mdict={'A': A,'nodes':obj_arr})
        nx.write_gexf(G_temp, 'lpus/v1_2_lpu_%s.gexf' % lpu)

        

if __name__ == "__main__":
    main()





