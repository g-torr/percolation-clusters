import numpy as np
import networkx as nx
import networkx.algorithms.bipartite as bipartite

import argparse
import sys
from numba import njit,jit
from numba.typed import Dict,List
import numba as nb
from scipy.special import factorial
sys.path.append('lib/')
from degree_seq_bipartite import bipartite_degree_seq
import os
import time
import pickle
from scipy.stats import pareto
def save_obj(obj, name ):
    with open('dic-'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def Convert_typed(tup):
    ''' for undirected edges. I construct the neighbouring relationship.
    Returns:
        di: dictionary.
    Notes:
    For link (a,b)  it returns the di[a]+=[b], and di[b]+=[a]
    '''
    d_typed = Dict.empty(
        key_type=nb.types.int64,
        value_type=nb.types.int64[:],
    )
    for a, b in tup:
        if a not in d_typed:
            d_typed[a] = np.array([b])
        else:
            d_typed[a] = np.append(d_typed[a], b)
        if b not in d_typed:
            d_typed[b] = np.array([a])
        else:
            d_typed[b] = np.append(d_typed[b], a)

    '''     
    for key, value in di.items():
        d_typed[key] = np.array(value,dtype=int) 
    '''
    return d_typed


def dict_to_typed_states(neigh):
    ''' for undirected edges. I construct the neighbouring relationship.
    Returns:
        di: dictionary.
    Notes:
    For link (a,b)  it returns the di[a]+=[b], and di[b]+=[a]
    '''
    d_typed = Dict.empty(
        key_type=nb.types.int64,
        value_type=nb.types.float32[:],
    )

    for key, value in neigh.items():
        d_typed[key] = np.ones(len(value), dtype=np.float32)
    return d_typed


def Convert(tup, di=None):
    ''' for undirected edges. I construct the neighbouring relationship.
    Returns:
        di: dictionary.
    Notes:
    For link (a,b)  it returns the di[a]+=[b], and di[b]+=[a]
    '''
    if di == None:
        di = {}
    for a, b in tup:
        di.setdefault(a, []).append(b)
        di.setdefault(b, []).append(a)
    di = {key: np.array(value, dtype=int) for key, value in di.items()}
    return di

def preparation(BG, usenumba=True):
    if usenumba:
        neigh = Convert_typed(list(BG.edges()))
        states = dict_to_typed_states(neigh)
    else:
        neigh = Convert(BG.edges())
        states = {key:np.ones(len(value)) for key,value in neigh.items()}#states[x][y] indicates the state of node "x cavity neigh[x][y]""
    return neigh,states
#states = {key:np.ones(len(value)) for key,value in neigh.items()}#states[x][y] indicates the state of node "x cavity neigh[x][y]""
@njit()# parallel version is way slower unfortunately
def single_instance_numba(neigh,states,N_a,N_b,p,N_iterations=20):
    individuals = set(range(N_a)).intersection(set(neigh.keys()))
    clusters = set(range(N_a,N_a+N_b)).intersection(set(neigh.keys()))
    #start = time.time()
    for count in range(N_iterations):    #solving cavity equations through forward dynamics
        err = 0
        for j in individuals:#solving 3rd equation
                for idx,mu in enumerate(neigh[j]):
                    new=1.0
                    for nu in set(neigh[j])-{mu}:
                        new*=1-states[nu][neigh[nu]==j][0]
                    err = max(abs(states[j][idx]-1+new),err)
                    states[j][idx]= 1-new
        for mu in clusters: #running 2nd equation
                for idx,i in enumerate(neigh[mu]):
                    new=1.
                    for j in set(neigh[mu])-{i}:
                        new*=(1-p*states[j][neigh[j]==mu][0])
                    err = max(abs(1-new-states[mu][idx]),err)
                    states[mu][idx]= 1-new
        if err<0.001:
            print('Exit after ',count,' iterations')
            break
    if count==N_iterations-1:
        print('Attention! Maximum number of iteration reached')

    #end = time.time()
    #print('it took',end-start,' seconds')

    risk ={}
    for i in individuals:#solving 3rd equation
        new=1
        for mu in set(neigh[i]):
            if len(states[mu][neigh[mu]==i])!=1:#sanity check
                print('mu=',mu,'i=',i)
                raise ValueError('Sanity check failed!')
            new*=(1-states[mu][neigh[mu]==i][0])
        risk[i]= 1-new
    return risk

def single_instance(neigh,states,N_a,N_b,p):
    individuals = set(range(N_a)).intersection(set(neigh.keys()))
    clusters = set(range(N_a,N_a+N_b)).intersection(set(neigh.keys()))
    #start = time.time()
    for count in range(N_iterations):    #solving cavity equations through forward dynamics
        for j in individuals:#solving 3rd equation
                for idx,mu in enumerate(neigh[j]):
                    new=1.
                    for nu in set(neigh[j])-{mu}:
                        new*=(1-states[nu][neigh[nu]==j])
                    states[j][idx]= 1-new
        for mu in clusters: #running 2nd equation
                for idx,i in enumerate(neigh[mu]):
                    new=1.
                    for j in set(neigh[mu])-{i}:
                        new*=(1-p*states[j][neigh[j]==mu])
                    states[mu][idx]= 1-new
    #end = time.time()
    #print('it took',end-start,' seconds')

    risk ={}
    for i in individuals:#solving 3rd equation
        new=1
        for mu in set(neigh[i]):
            new*=(1-states[mu][neigh[mu]==i])
        risk[i]= 1-new
    return risk

def save_obj(obj):
    directory = 'data'
    if not os.path.exists(directory):
        os.makedirs(directory)
    name = "unstructured.pkl"
    with open(directory + '/dic-' + name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def main():
    parser = argparse.ArgumentParser(
        description='Single instance cavity for unstructured contagion in undirected bipartite network.  Shifted Poisson connectivity. Default run with numba for maximum performance. \n'
                    'However, numba is under development so it can give problems with future updates. Run with the --no-rumba option in case')
    parser.add_argument("--no_numba", help="Avoid use of the numba package",  action="store_true")
    args = parser.parse_args()
    usenumba = not args.no_numba
    N_a = 100_000# n people
    N_b = 100_000 #n clusters
    a_mean = 4. #average degree of the distribution
    b_mean = (a_mean )* N_a / N_b
    p = 0.1 #probability to transmit a disease
    #----END of parameter definition-------
    #if you want a different distribution, replace shifted_poisson with any random generator in np.random
    aseq,bseq = bipartite_degree_seq(N_a,N_b,'shifted_poisson','shifted_poisson',{'lam':a_mean},{'lam':b_mean})
    #aseq,bseq = bipartite_degree_seq(N_a,N_b,'poisson','poisson',{'lam':a_mean},{'lam':N_a*a_mean/N_b})
    BG = bipartite.generators.configuration_model(aseq,bseq)
    BG = nx.Graph(BG)#convert multilinks to simple
    neigh, states = preparation(BG, usenumba)
    start = time.time()
    if usenumba:
        risk = single_instance_numba(neigh, states,N_a,N_b,p)
    else:
        risk = single_instance(neigh, states,N_a,N_b,p)
    stop = time.time()
    print(' It tooks ',stop-start,' s to execute.')
    dic ={"N_a":N_a,"N_b":N_b,'aseq':aseq,'bseq':bseq,'risk':dict(risk)}
    save_obj(dic)
if __name__ == '__main__':
    main()
