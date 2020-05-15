import math
import random
import networkx 
from functools import reduce
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys
version=sys.version_info[:2]
import datetime
print("module loaded at ",datetime.datetime.now())

def choices(population,k):
		'''Backward compatibility for python 3.5'''
		total = len(population)
		return [population[int(random.random() * total)] for i in range(k)] 
def shifted_poisson(lam = 1,size = None):
    '''
    Draw samples from a 1+Poisson distribution.

    The Poisson distribution is the limit of the binomial distribution
    for large N.

    Parameters
    ----------
    lam : float or array_like of floats
        Expectation of interval, must be >= 1. A sequence of expectation
        intervals must be broadcastable over the requested size.
        Poisson distribution is extracted with average degree lam-1
    size : int or tuple of ints, optional
        Output shape.  If the given shape is, e.g., ``(m, n, k)``, then
        ``m * n * k`` samples are drawn.  If size is ``None`` (default),
        a single value is returned if ``lam`` is a scalar. Otherwise,
        ``np.array(lam).size`` samples are drawn.

    Returns
    -------
    out : ndarray or scalar
        Drawn samples from the parameterized Poisson distribution.
    '''

    
    return np.random.poisson(lam = lam-1,size = size)+1
np.random.shifted_poisson = shifted_poisson
if version[0]<3:
	raise Error("Python 2 not supported")
if version[1]<6:
	print("ciao")
	random.choices=choices

		
		
def bipartite_degree_seq(N_a,N_b,pmf_a,pmf_b,a_params,b_params):
    '''
    generate the degree sequences according to the degree distribution specified.
    Accepted degree distributions are list at np.random.__all__
    Parameters:
    N_a: int
        number of nodes in layer a
    N_b: int
        number of nodes in layer b        
    pmf_a: string
        name of the distribution from which degree sequence is sampled, for nodes a
    pmf_b: string
        name of the distribution from which degree sequence is sampled, for nodes b
    a_params: dict
        parameters of the distribution pmf_a
    b_params: dict
        parameters of the distribution pmf_b
    Returns:
    aseq: np.array
        degree sequence layer a
    bseq: np.array
        degree sequence layer b
    Usage:
    N_a = 10000
    N_b = 1000
    a_mean = 1
    a,b = bipartite_degree_seq(N_a,N_b,'poisson','poisson',{'lam':a_mean},{'lam':N_a*a_mean/N_b})

    '''
    print(a_params)
    afunc = getattr(np.random,pmf_a)
    aseq =afunc(**a_params,size = N_a) 
    bseq =afunc(**b_params,size = N_b) 
    diff =np.sum(aseq)- np.sum(bseq) 
    if abs(diff/np.sum(aseq))>0.01:
        raise ValueError('It is likely that statistically the sum of the two degree sequence differs')
    while diff != 0:
        n_new_samples = len(aseq) #int(abs(diff) /a_mean)
        new_samples =afunc(**a_params,size = n_new_samples)
        for proposed in new_samples:
            i = np.random.choice(len(aseq))
            k = aseq[i]
            if abs(diff -k+ proposed )<abs(diff):
                aseq[i] = proposed
                diff = diff -k+ proposed
    return aseq,bseq


