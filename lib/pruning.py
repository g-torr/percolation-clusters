import numpy as np
import scipy
import sys
import itertools
from multiprocessing import Pool
import pickle
import random
import networkx as nx
sys.path.insert(0, "./lib")# add the library folder to the path I look for modules
import component as co

def out_nonlinear_remove(BG, Genes, TFs, genes_rm):
    '''Return the subgraph in the SCC, and then removes TFs whose in degree is not the same as initially '''
    BG_filtered = BG.copy()
    BG_tf = dict(BG.in_degree(TFs))  # dictionary of in degree for TF in the bipartite graph
    BG_filtered.remove_nodes_from(genes_rm)
    count = 0
    while True:
        Gc_nodes = co.out_component(BG_filtered)
        Gc = BG_filtered.subgraph(Gc_nodes).copy()
        if len(Gc_nodes) == len(BG_filtered):
            print("number of deletion iterations:", count)
            return BG_filtered, count
        count += 1

        Gc_G = set(Gc.nodes()).intersection(Genes)  # set of genes in the giant cluster
        Gc_tf = dict(Gc.in_degree(TFs))  # dictionary of in degree for TF in the SCC
        filtered = {v for v, k in Gc_tf.items() if
                    BG_tf[v] == k}  # filtered dictionary of SCC if TF has same indegree as initially
        if len(filtered) == 0:
            return nx.empty_graph(), count
        BG_filtered = BG_filtered.subgraph(set(filtered).union(Gc_G)).copy()



def filter_out(BG):
    '''It extract the out subgraph'''
    nodes_out = co.out_component(BG)
    Gc = BG.subgraph(nodes_out).copy()
    Genes = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 0]
    N1 = len(Genes)
    TFs = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 1]
    N2 = len(TFs)
    in_deg = scipy.sparse.csr_matrix(list(dict(Gc.in_degree(TFs)).values()))
    return Genes, TFs, in_deg

