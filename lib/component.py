'''This module computes in-,out- component and SCC. Contains functions:
in_component(G)
out_component(G)
SCC(G) 
---
There are 2  classes for numerical calulations:
-num_Poisson(d,c,alpha)
-num_corrected(d,c,alpha)
'''
from pylab import *
from scipy.optimize import fsolve
import networkx as nx
def maximum_strongly_connected_component_subgraph(BG):
    ''' Returns the subgraph of strongly connected component of BG'''
    SCC_nodes=max([c for c in nx.strongly_connected_components(BG)],key=len)
    return BG.subgraph(SCC_nodes)
def in_component(G):
    ''' Returns the set of nodes in the in-component of G'''
    Gc=maximum_strongly_connected_component_subgraph(G)
    i=list(Gc)[0]
    #i=next(Gc.nodes_iter()) Networkx 1.
    return nx.ancestors(G,i).union({i})
def out_component(G):
    ''' Returns the set of nodes in the out-component of G'''
    Gc=maximum_strongly_connected_component_subgraph(G)
    i=list(Gc)[0]
    #i=next(Gc.nodes_iter())
    return nx.descendants(G,i).union({i})
def SCC(G):
    ''' Returns the set of nodes in the strongly connected component of G'''
    Gc=maximum_strongly_connected_component_subgraph(G)
    return set(Gc.nodes())
class num_Poisson:
	def __init__(self,d,c,alpha):
		self.d=d
		self.c=c
		self.alpha=alpha

	def __get__in_component(self):
		d=self.d
		c=self.c
		alpha=self.alpha
		gene=lambda g,d,c: 1-exp(d*(-1+exp(-c * g)))-g
		TF=lambda  g,d,c: 1-exp(c*(-1+exp(-d*g)))- g
		return array([fsolve(gene,1,(d,c)),fsolve(TF,1,(d,c))])
	def __get__out_component(self):
		def gene(g,d,c,alpha):
			return 1-exp(c*(-1+exp(-d * g/alpha))*alpha)-g     
		def TF(g,d,c,alpha):
			return 1-exp(d*(-1+exp(-c * g*alpha))/alpha)-g
		d=self.d
		c=self.c
		alpha=self.alpha

		return array([fsolve(gene,1,(d,c,alpha)),fsolve(TF,1,(d,c,alpha))])  
	def __get__SCC(self):
		d=self.d
		c=self.c
		alpha=self.alpha
		return self.in_component*self.out_component
	in_component = property(__get__in_component)
	out_component = property(__get__out_component)
	SCC = property(__get__SCC)
		
class num_corrected:
	def __init__(self,d,c,alpha):
		self.d=d
		self.c=c
		self.alpha=alpha

	def __get__in_component(self):
		d=self.d
		c=self.c
		alpha=self.alpha

		def gene(g,d,c):
			return 1-exp(d*(-1+exp(-c*g))-c*g) -g       
		def TF(g,d,c):
			return 1-exp(-(c)*exp(-d*g)*(-1+exp(d*g)+g))-g
		d_out=d+where(alpha>1,alpha-1,0)
		c_out=c+1/alpha
		return array([fsolve(gene,1,(d_out,c_out)),fsolve(TF,1,(d_out,c_out))])
	def __get__out_component(self):
		d=self.d
		c=self.c
		alpha=self.alpha
		def gene(g,d,c):
			return exp(-c - d*g - c*(-1 + g)*exp(-(d*g)))*(-1 + g + exp(c + d*g + c*(-1 + g)*exp(-(d*g))))-g     
		def TF(g,d,c):
			return exp(-d - c*g - d*(-1 + g)*exp(-(c*g)))*(-1 + g + exp(d + c*g + d*(-1 + g)*exp(-(c*g))))-g
		d_in=c*alpha
		c_in=d/alpha+where(alpha<1,1/alpha-1,0)
		return array([fsolve(gene,1,(c_in,d_in)),fsolve(TF,1,(c_in,d_in))])# I exchange c and d on purpose. In Mathematica it stands for the alpha=1 case, so c is the indegree of TF, and d is the indegree of Genes
	def __get__SCC(self):
		d=self.d
		c=self.c
		alpha=self.alpha
		return self.in_component*self.out_component
	in_component = property(__get__in_component)
	out_component = property(__get__out_component)
	SCC = property(__get__SCC)
            
