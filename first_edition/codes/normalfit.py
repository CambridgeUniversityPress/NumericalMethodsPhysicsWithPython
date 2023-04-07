# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from gauelim_pivot import gauelim_pivot
import numpy as np

def generatedata(N):
    np.random.seed(45379)
    dataxs = np.linspace(0,9,N)
    datays = 2 + 5*np.sin(dataxs) + 0.3*np.random.randn(N)
    datasigs = 0.2*np.abs(np.random.randn(N))
    return dataxs, datays, datasigs

def phi(n,k,x):
    if n==5:
        val = x**k
    elif n==2:
        val = 1. if k==0 else np.sin(x)
    return val

def normalfit(dataxs,datays,datasigs,n):
    N = dataxs.size
    A = np.zeros((N,n))
    for k in range(n):
        A[:,k] = phi(n,k,dataxs)/datasigs
    bs = datays/datasigs

    cs = gauelim_pivot(A.T@A, A.T@bs)
    chisq = np.sum((bs - A@cs)**2)
    return cs, chisq

if __name__ == '__main__':
    dataxs, datays, datasigs = generatedata(8)
    for n in (5, 2):
        cs, chisq = normalfit(dataxs, datays, datasigs, n)
        print(cs)
        print(chisq/(dataxs.size-cs.size)); print("")
