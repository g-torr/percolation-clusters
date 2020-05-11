import argparse
from directed_random_graph import *# this is the function I wrote to create a random bipartite
def create_network(args):
	''' This function creates the network from input values. Makes sanity checks on the input parameters'''
	if (len(args.out_degree)>2)or(len(args.out_degree)==1) :
		raise ValueError(" Out degree parameter has to be  two positive float, which corresponds to c and d. Open help page for more info ")
	
	c,d=args.out_degree
	if (c<0 or d<0):
		raise ValueError(" Out degree parameter has to be  two positive float, which corresponds to c and d. Open help page for more info ")

	N1=args.N1 # number genes
	N2=args.N2 # number TFs
	corrected=not args.simple_Poisson
	print("creating the network with ",N1," genes. ",N2," TF. Corrected=", corrected, ". c= ",c,". d= ",d )
	BG, C =create_bipartite(N1,N2,d,c)
	print("network completed.")
	if corrected:
		return C,c,d,N1,N2
	else:
		return BG,c,d,N1,N2 
