import math
import random
import networkx 
from functools import reduce
import networkx as nx
import numpy.random
import matplotlib.pyplot as plt
import sys
version=sys.version_info[:2]
import datetime
print("module loaded at ",datetime.datetime.now())
def choices(population,k):
		'''Backward compatibility for python 3.5'''
		total = len(population)
		return [population[int(random.random() * total)] for i in range(k)] 

if version[0]<3:
	raise Error("Python 2 not supported")
if version[1]<6:
	print("ciao")
	random.choices=choices

		
		
 
def create_bipartite(N1,N2,dout,cout):
    ''' e' la funzione che ho scritto io. This is a faster version that reduce the numer of random permuations evaluated '''

    N=min(N1,N2)

    BG = random_graph(N1, N2, (dout/N2,cout/N1),directed=True)
    top_nodes = [n for n,d in BG.nodes(data=True) if d['bipartite']==0]
    bottom_nodes = [n for n,d in BG.nodes(data=True) if d['bipartite']==1]
    F=BG.copy()
    new_edges=list(zip(numpy.random.permutation(top_nodes[:N]),bottom_nodes[:N]))# add gene to TF links
    if N==N2:# i.e. if genes are more than TF
        new_edges+=list(zip(top_nodes[N:],random.choices(bottom_nodes,k=abs(N1-N2))))
    else:
        new_edges+=list(zip(random.choices(top_nodes,k=abs(N1-N2)),bottom_nodes[N:]))        
    
    F.add_edges_from(new_edges)
    
    new_edges=list(zip(random.choices(bottom_nodes,k=N1),top_nodes))# add TF to gene links. From a random sampling of TF, one take the  Genes to connect to 
    F.add_edges_from(new_edges)
    return BG,F
def random_graph(n, m, prob, seed=None, directed=False):
    """L'ho scritta io. Nodi possono avere degree 0
    Return a bipartite random graph.

    This is a bipartite version of the binomial (Erdos-Renyi) graph.

    Parameters
    ----------
    n : int
        The number of nodes in the first bipartite set.
    m : int
        The number of nodes in the second bipartite set.
    prob : float (undirected graph), 2-dim tuple (for directed graph) 
        Probability for edge creation.
        p[0] is the probability of set 1 to set 2
        p[1] is the probability of set 2 to set 1
    seed : int, optional
        Seed for random number generator (default=None). 
    directed : bool, optional (default=False)
        If True return a directed graph 
      
    Notes
    -----
    This function is not imported in the main namespace.
    To use it you have to explicitly import the bipartite package.

    The bipartite random graph algorithm chooses each of the n*m (undirected) 
    or 2*nm (directed) possible edges with probability p.

    This algorithm is O(n+m) where m is the expected number of edges.
    
    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.

    See Also
    --------
    gnp_random_graph, configuration_model

    References
    ----------
    .. [1] Vladimir Batagelj and Ulrik Brandes, 
       "Efficient generation of large random networks",
       Phys. Rev. E, 71, 036113, 2005.
    """
    if (type(prob)==tuple) and (directed):
        if len(prob)!=2:
             raise ValueError("parameter prob should be a 2-dim tuple")
    else:
        raise ValueError(" if you set a directed graph, you need a tuple of probabilites ( they may be the same, e.g. prob=(0.1,0.1))")
    
    G=nx.Graph()
    G=_add_nodes_with_bipartite_label(G,n,m)
    if directed:
        G=nx.DiGraph(G)
    G.name="fast_gnp_random_graph(%s,%s,%s)"%(n,m,prob)

    if not seed is None:
        random.seed(seed)

    if directed:
        p=prob[0]
    else:
        p=prob
        if p <= 0:
            return G
        if p >= 1:
            return nx.complete_bipartite_graph(n,m)

    lp = math.log(1.0 - p)  

    v = 0 
    w = -1
    while v < n:
        lr = math.log(1.0 - random.random())
        w = w + 1 + int(lr/lp)
        while w >= m and v < n:
            w = w - m
            v = v + 1
        if v < n:
            G.add_edge(v, n+w)

    if directed:
        p1=prob[1]
        # use the same algorithm to 
        # add edges from the "m" to "n" set
        lp1 = math.log(1.0 - p1)
        v = 0 
        w = -1
        while v < n:
            lr = math.log(1.0 - random.random())
            w = w + 1 + int(lr/lp1)
            while  w>= m and v < n:
                w = w - m
                v = v + 1
            if v < n:
                G.add_edge(n+w, v)

    return G
def _add_nodes_with_bipartite_label(G, lena, lenb):
    G.add_nodes_from(range(0,lena+lenb))
    b=dict(zip(range(0,lena),[0]*lena))
    b.update(dict(zip(range(lena,lena+lenb),[1]*lenb)))
    if networkx.__version__>'2':
         nx.set_node_attributes(G,b,'bipartite')
    else:
         nx.set_node_attributes(G,'bipartite',b)         
    return G

def plot_bipartite(BG,nodelist=None,nodecolor="r"):
    '''nodelist: (optional) gives the set of nodes to plot''' 
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
def maximum_strongly_connected_component_subgraph(BG):
    SCC_nodes=max([c for c in nx.strongly_connected_components(BG)],key=len)
    return BG.subgraph(SCC_nodes)

