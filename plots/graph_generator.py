import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools

def draw_and_plot(G,node_colour="blue",edge_colour="black",with_labels=False,node_size=300):
    #https://networkx.github.io/documentation/networkx-1.9/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html
    nx.draw_networkx(G,node_color=node_colour,edge_color=edge_colour,with_labels=with_labels,node_size=node_size)
    plt.draw()
    plt.show()

# randomly partitioned graph 
# random_partition_graph(sizes, p_in, p_out, seed=None, directed=False)
node_size = 1
number_clusters = 8
nodes = [node_size for _ in range(number_clusters)]
G = nx.random_partition_graph(nodes,0.2,1.00,directed=False)


# G = nx.powerlaw_cluster_graph(200,5,0.0)
G = nx.relaxed_caveman_graph(8, 16, 0.05)

#draw_and_plot(G)

# custom graph (hyperedges)
#G = nx.DiGraph()
#nodes = [x for x in range(16)]
#G.add_nodes_from(nodes)
#edges0 = [(15,0),(15,8)]
#G.add_edges_from(edges0,color='black')
#edges1 = [(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,7)]
#G.add_edges_from(edges1,color='r')
#edges2 = [(8,9),(8,10),(8,11),(8,12),(8,13),(8,14)]
#G.add_edges_from(edges2,color='g')
#edge_colours = [G[u][v]['color'] for u,v in G.edges()]
#draw_and_plot(G,with_labels=True,node_colour='white',edge_colour=edge_colours,node_size=700)

# custom graph (load balance)
#G = nx.DiGraph()
#nodes = [x for x in range(8)]
#G.add_nodes_from(nodes)
#pos = nx.spring_layout(G)
#edges0 = [(1,0),(2,0),(3,0),(4,0),(5,1),(5,2),(5,3),(5,4),(6,0),(6,4),(6,5),(7,0),(7,5),(7,1)]
#G.add_edges_from(edges0,color='black')
#edge_colours = [G[u][v]['color'] for u,v in G.edges()]
#labels={}
#labels[0]='7'
#labels[1]='3'
#labels[2]='2'
#labels[3]='2'
#labels[4]='3'
#labels[5]='3'
#labels[6]='1'
#labels[7]='1'
#nx.draw_networkx(G,pos=pos,node_color='white',edge_color=edge_colours,with_labels=False,node_size=2000)
#nx.draw_networkx_labels(G,pos,labels,font_size=16)
#plt.draw()
#plt.show()

# custom graph (synfire)
G = nx.Graph()
nodes1 = [x for x in range(0,10)]
nodes2 = [x for x in range(10,20)]
nodes3 = [x for x in range(20,30)]
nodes4 = [x for x in range(30,40)]
G.add_nodes_from(nodes1)
G.add_nodes_from(nodes2)
G.add_nodes_from(nodes3)
G.add_nodes_from(nodes4)
edges1 = list(itertools.product(nodes1, nodes2))
edges2 = list(itertools.product(nodes2, nodes3))
edges3 = list(itertools.product(nodes3, nodes4))
G.add_edges_from(edges1,color='black')
G.add_edges_from(edges2,color='black')
G.add_edges_from(edges3,color='black')
edge_colours = [G[u][v]['color'] for u,v in G.edges()]
draw_and_plot(G,with_labels=False,node_colour='blue',edge_colour=edge_colours,node_size=700)