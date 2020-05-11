import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.bipartite as bipartite
import random
#generators.configuration_model as configuration_model
def _add_nodes_with_bipartite_label(G, lena, lenb):
    G.add_nodes_from(range(0,lena+lenb))
    b=dict(zip(range(0,lena),[0]*lena))
    b.update(dict(zip(range(lena,lena+lenb),[1]*lenb)))
    if nx.__version__>'2':
         nx.set_node_attributes(G,b,'bipartite')
    else:
         nx.set_node_attributes(G,'bipartite',b)         
    return G

def pair_stubs(nodes_from,nodes_to,kseq):
    '''Returns the list of edges.
    Parameters
    ----------
    nodes_from : list
        List of nodes from which links depart.
    nodes_to : list
        List of nodes where links arrive.
    kseq : list
        List of out-degree sequence of "nodes_from"


    It creates the stubs "nodes_from"-> "nodes_to".
    "nodes_from" are then paired to "nodes_to". For each node in node_from,
    its neighbours in set "nodes_to" are randomly sampled without replacement. '''
    # build lists of degree-repeated vertex numbers in the s
    stubs = []
    m = min(nodes_from) # m=0 if "nodes_from" is set A, m = len(aseq) if n "nodes_from" is set B
    stubs.extend([[v] * kseq[v-m] for v in nodes_from])
    stubs_from = []
    stubs_from = [x for subseq in stubs for x in subseq]
    # build list of arrival vertex numbers in b
    stubs_to = [random.sample(nodes_to,k_out) for k_out in kseq]
    stubs_to = [x for subseq in stubs_to for x in subseq]
    return stubs_from,stubs_to

def configuration_model(aseq, bseq, create_using=None, seed=None):
    """Returns a random bipartite graph from two given degree sequences.
    Directed graph with aseq and bseq being the out-degree.

    Parameters
    ----------
    aseq : list
       Out-degree sequence for node set A.
    bseq : list
       Out-degree sequence for node set B.
    create_using : NetworkX graph instance, optional
       Return graph of this type.
    seed : integer, random_state, or None (default)
        Indicator of random number generation state.
        See :ref:`Randomness<randomness>`.

    The graph is composed of two partitions. Set A has nodes 0 to
    (len(aseq) - 1) and set B has nodes len(aseq) to (len(bseq) - 1).
    Nodes from set A are connected to nodes in set B by choosing
    randomly from the possible free stubs, one in A and one in B.

    Notes
    -----
    The sum of the two sequences must be equal: sum(aseq)=sum(bseq)
    If no graph type is specified use MultiGraph with parallel edges.
    If you want a graph with no parallel edges use create_using=Graph()
    but then the resulting degree sequences might not be exact.

    The nodes are assigned the attribute 'bipartite' with the value 0 or 1
    to indicate which bipartite set the node belongs to.

    This function is not imported in the main namespace.
    To use it use nx.bipartite.configuration_model
    """
    G = nx.empty_graph(0, create_using, default=nx.MultiGraph)

    # length and sum of each sequence
    lena = len(aseq)
    lenb = len(bseq)
    suma = sum(aseq)
    sumb = sum(bseq)
    nodes_a = range(lena)
    nodes_b = range(lena, lenb + lena)

    G = _add_nodes_with_bipartite_label(G, lena, lenb)

    if len(aseq) == 0 or max(aseq) == 0:
        return G  # done if no edges
    # a->b links
    astubs, bstubs = pair_stubs(nodes_a, nodes_b, aseq)
    G.add_edges_from(zip(astubs, bstubs))
    # b->a links
    bstubs, astubs = pair_stubs(nodes_b, nodes_a, bseq)
    G.add_edges_from(zip(bstubs, astubs))

    G.name = "bipartite_configuration_model"
    return G


