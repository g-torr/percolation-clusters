import networkx as nx
import itertools
import random
import numpy as np
import numpy.random as npr

def _add_nodes_with_bipartite_label(G, lena, lenb):
    G.add_nodes_from(range(0,lena+lenb))
    b=dict(zip(range(0,lena),[0]*lena))
    b.update(dict(zip(range(lena,lena+lenb),[1]*lenb)))
    if nx.__version__>'2':
         nx.set_node_attributes(G,b,'bipartite')
    else:
         nx.set_node_attributes(G,'bipartite',b)
    return G
def pair_stubs_regulatory(nodes_from, nodes_to, kseq):
    '''Returns the list of edges.
    Parameters
    ----------
    nodes_from : list
        List of nodes from which links depart.
    nodes_to : list
        List of nodes where links arrive.
    kseq : np.array
        Array of out-degree sequence of "nodes_from"

    Returns:
    from: np.array
        Array of node ids sampled from "nodes_from"
    to: np.array
        Array of nodes ids sampled from "nodes_to"


    It creates the stubs "nodes_from"-> "nodes_to".
    "nodes_from" are then paired to "nodes_to". Links are
    generated such that  in-degree of "nodes_to" is >=1.
    ----NOTE---
    It is possible to have multi-links.

    Stubs are generated as following:
    - reg stubs: exactly len(nodes_to) stubs are sampled among all possible
    - nonreg stubs: all the other stubs
    nonreg stubs never are paired to the same nodes,therefore multi links are excluded in this set.
    However, there is no check that reg and nonreg stubs repeat a link.
    '''

    # build lists of degree-repeated vertex numbers in the s
    if sum(kseq) < len(nodes_to):
        raise ValueError('Degree sequence is not compatible with regulatory couplings')
    r, r_seq = reg_stubs(nodes_from, nodes_to, kseq)
    nodes_reg_from = [[x] * x_seq for x, x_seq in zip(r, r_seq)]  # to impose the in-degree is at least 1
    nodes_reg_from = list(itertools.chain.from_iterable(nodes_reg_from))
    nodes_reg_to = npr.permutation(nodes_to)  # to impose the in-degree is at least 1
    # --- now deals with non-reg stubs
    stubs = []
    m = min(nodes_from)  # m=0 if "nodes_from" is set A, m = len(aseq) if n "nodes_from" is set B

    stubs_from, stubs_kseq = non_reg_stubs(nodes_from, nodes_to, kseq, r, r_seq)
    # return stubs_from,stubs_kseq,r, r_seq
    stubs.extend([[v] * k for v, k in zip(stubs_from, stubs_kseq)])
    stubs_from = [x for subseq in stubs for x in subseq]
    # build list of arrival vertex numbers in b
    stubs_to = [random.sample(nodes_to, k_out) for k_out in stubs_kseq]
    stubs_to = [x for subseq in stubs_to for x in subseq]
    return np.append(nodes_reg_from, stubs_from), np.append(nodes_reg_to, stubs_to)




def reg_stubs(nodes_from, nodes_to, kseq):
    '''
    We sample the #stubs= len(nodes_to) from the set "nodes_from".
    We use this set of nodes and degrees to create a regulatory network,
    such that the in-degree of nodes_to is >=1".
        Parameters
    ----------
    nodes_from : list
        List of nodes from which links depart.
    nodes_to : list
        List of nodes where links arrive.
    kseq : np.array
        Array of out-degree sequence of "nodes_from"

    Returns:
    r: list.
        List of node ids sampled from "nodes_from"
    r_seq: np.array
        Array of degrees of nodes in R. It holds sum(f_seq)== len(nodes_to)
    They are the nodes whose stubs are used for regulatory links.
    ----
    NOTE:
    Multilinks are possible.'''

    m = min(nodes_from)
    a = npr.choice(list(itertools.chain.from_iterable([[v] * kseq[v - m] for v in nodes_from])), size=len(nodes_to),
                   replace=False)
    r, r_seq = np.unique(a, return_counts=True)
    return r, r_seq


def non_reg_stubs(nodes_from, nodes_to, kseq, r, r_seq):
    '''Opposite job to "reg_stubs" function. It creates the list of nodes in "nodes_from"
    whose stubs are not (fully) used for the regulatory part.
        Parameters
    ----------
    nodes_from : list
        List of nodes from which links depart.
    nodes_to : list
        List of nodes where links arrive.
    kseq : np.array
        Array of out-degree sequence of "nodes_from"

    Returns:
        nr: list.
            List of nodes ids selected from "nodes_from".
        nr_seq: np.array()
            Array of stubs departing from corresponding nodes.
        '''
    m = min(nodes_from)  # m=0 if "nodes_from" is set A, m = len(aseq) if n "nodes_from" is set B
    diff = kseq[np.array(
        r) - m] - r_seq  # ==0 means that the stubs of that node have been used entirely for regulatory stubs
    nr = list(np.append(np.array(r)[diff > 0], np.array(nodes_from)[~np.isin(nodes_from,
                                                                             r)]))  # concatenate the list of nodes partially used for regulation, with the nodes that have not been used for regulation at all
    nr_seq = np.append(diff[diff > 0], kseq[~np.isin(nodes_from, r)])
    return nr, nr_seq

def configuration_model(aseq, bseq, create_using=None, seed=None):
    """Returns a random bipartite graph from two given degree sequences.

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
    In degree of nodes is 1+k where k is a random variable from Poisson
    distribution.

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
    G = nx.empty_graph(0, create_using, default=nx.MultiDiGraph)

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
    astubs, bstubs = pair_stubs_regulatory(nodes_a, nodes_b, aseq)
    G.add_edges_from(zip(astubs, bstubs))
    # b->a links
    bstubs, astubs = pair_stubs_regulatory(nodes_b, nodes_a, bseq)
    G.add_edges_from(zip(bstubs, astubs))

    G.name = "bipartite_configuration_model"
    return G

def mirror_stubs(nodes_from,nodes_to,kseq):
    '''
    Create in- stubs  among the set "nodes_to"
    '''
    m = min(nodes_from)  # m=0 if "nodes_from" is set A, m = len(aseq) if n "nodes_from" is set B
    stubs_from = list(itertools.chain.from_iterable([[v] * kseq[v - m] for v in nodes_from]))
    bseq_in = kseq.copy()
    m = min(nodes_to)
    np.random.shuffle(bseq_in)#in-degree sequence of "nodes_to"
    stubs_to = list(itertools.chain.from_iterable([[v] * bseq_in[v - m] for v in nodes_to]))
    np.random.shuffle(stubs_from)
    return  stubs_from,stubs_to

def configuration_model_mirroring(aseq, bseq, create_using=None, seed=None):
    """Returns a random bipartite graph from two given degree sequences.

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
    The in-degree sequence of nodes in B is the same but shuffled of
    the out-degree in A.
    Nodes from set B are connected to nodes in set A by choosing
    randomly from the possible free stubs, one in B and one in A.
    In degree of nodes in A is 1+k where k is a random variable from Poisson
    distribution.

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
    G = nx.empty_graph(0, create_using, default=nx.MultiDiGraph)

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
    if lena != lenb:
        raise ValueError('This version requires the same number of nodes on the 2 layers. Hope to improve')
    # a->b links
    astubs, bstubs = mirror_stubs(nodes_a, nodes_b, aseq)
    G.add_edges_from(zip(astubs, bstubs))
    # b->a links
    bstubs, astubs = pair_stubs_regulatory(nodes_b, nodes_a, bseq)
    G.add_edges_from(zip(bstubs, astubs))

    G.name = "bipartite_configuration_model"
    return G
