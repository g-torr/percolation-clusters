import networkx as nx
from pylab import *
import matplotlib.pyplot as plt
def plot_bipartite(BG,nodelist=None,nodecolor="r"):
    print("Edges from left to right are colored green. Right to left: blue, bidirectional: magenta.")
    if nodelist != None:
        BG=nx.subgraph(BG,nodelist)
    '''
    Plot bipartite graph, set 1 in x=0, set 2 in x=1.  
    Color for links from set 1 to set 2 in green,
    from set 2 to set 1 in blue,
    magenta if there is a double link between pair of nodes'''
    genes = [n for n,d in BG.nodes(data=True) if d['bipartite']==0]
    TFs = [n for n,d in BG.nodes(data=True) if d['bipartite']==1]
    pos = dict()
    pos.update( (n, (1, i)) for i, n in enumerate(genes) ) # put bottom_nodes  x=1
    pos.update( (n, (2, i)) for i, n in enumerate(TFs) ) # put top_nodes at x=2
    edgelist,edge_color=color_edges(BG)
    nx.draw_networkx(BG, pos=pos,with_labels=True,edgelist=edgelist,edge_color=edge_color,nodelist=nodelist,node_color=nodecolor)
    plt.text(1,len(genes),"G")
    plt.text(2,len(TFs),"TF")

def color_edges(BG):
    '''Color edges of bipartite. Nodes in graph need to have a "bipartite" label.
    This coloring relies on the assumption that the 2 sets of nodes are separated.
    Color for links from set 1 to set 2 in green, from set 2 to set 1 in blue, magenta if there is a double link between pair of nodes
    Input:
        BG: nx.Digraph with bipartite label in each node
    Returns:
        links: list . list of touples containing the node ids for a link
        colors: list . List of colors, same dimension of links '''
    
    node1=[]
    node2=[]
    for n,d in BG.nodes(data=True):
        if (len(d)==0 or list(d.keys())[0]!="bipartite"):
            raise TypeError("Graph is not a bipartite")
        if d['bipartite']==0:
            node1+=[n]
        elif d['bipartite']==1:
            node2+=[n]
        else:
            raise TypeError("bipartite label should contains only value 0 or 1, not "+str(d["bipartite"]))
    if max(node1)>min(node2):
        if (min(node1)<max(node2)):
            raise TypeError("nodes id of two classes of nodes should be separate, this is not the case here")
    colors=[]
    links=[l for l in BG.edges()]
    for a,b in BG.edges():
        if (b,a) in links:#if there are double links 
            colors+=["m"]
        elif a>b:
            colors+=["b"]
        elif b>a:
            colors+=["g"]
        else:
            raise ValueError("Suspected self node link, which should not happen in bipartite")
    return links,colors
