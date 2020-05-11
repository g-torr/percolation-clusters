from pylab import *
from scipy.stats import poisson
def Poissoninan_out(d,c,M):
    P=[]
    for k in  arange(max(data)):
        s=0
        for n in arange(10*M):
            s+=poisson.pmf(n, d)*poisson.pmf(k, n*c)
        P+=[s]
    return arange(1,M),P
    plot(P)
def Poissoninan_in(d,c,alpha,M):
    d_in=c*alpha
    c_in=d/alpha+where(alpha<1,1/alpha-1,0)
    P=[]
    for k in  arange(M):
        s=0
        for n in arange(10*M):
            s+=poisson.pmf(n, d_in)*poisson.pmf(k, n*c_in)
        P+=[s]
    return arange(1,M),P
    plot(P)
def corrected_Poissoninan_out(d,c,alpha,M):
    d_out=d+where(alpha>1,alpha-1,0)
    c_out=c+1/alpha
    P=[]
    for k in  arange(M):
        s=0
        for n in arange(10*M):# this sum should go run till infinity
            s+=poisson.pmf(n, d_out)*poisson.pmf(k, (n+1)*c_out)
        P+=[s]
    return arange(M),P
    plot(arange(max(data)),P,label="theory")
def corrected_Poissoninan_in(d,c,alpha,M):
    d_in=c*alpha
    c_in=d/alpha+where(alpha<1,1/alpha-1,0)
    P=[]
    for k in  arange(1,M):
        s=0
        for n in arange(k):
            s+=poisson.pmf(n, d_in)*poisson.pmf(k-n-1, (1+n)*c_in)
        P+=[s]
    return arange(1,M),P 
    plot( arange(1,max(data)),P,label="theory")

