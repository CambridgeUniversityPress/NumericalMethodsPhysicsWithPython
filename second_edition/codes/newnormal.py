# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from gauelim_pivot import gauelim_pivot
import numpy as np

def generatedata(N,a=0.,b=9,sts=(2,5,0.5,1)):
    sa, sb, sc, sd = sts
    np.random.seed(7921)
    data = np.zeros((3,N))
    data[0,:] = np.linspace(a,b,N)
    data[1,:] = sa + sb*np.sin(data[0,:])
    data[2,:] = sc + sd*np.random.random(N)
    data[1,:] += np.random.normal(0,data[2,:])
    return data

def phi(n,k,x):
    if n==5:
        val = x**k
    elif n==2:
        val = 1. if k==0 else np.sin(x)
    return val

def normalfit(data,phi,n):
    N = data.shape[1]
    A = np.zeros((N,n))
    for k in range(n):
        A[:,k] = phi(n,k,data[0,:])/data[2,:]
    bs = data[1,:]/data[2,:]

    cs = gauelim_pivot(A.T@A, A.T@bs)
    chisq = np.sum((bs - A@cs)**2)
    return cs, chisq

if __name__ == '__main__':
    data = generatedata(8)
    for n in (5, 2):
        cs, chisq = normalfit(data, phi, n)
        print(cs)
        print(chisq/(data.shape[1]-cs.size))
        print("")
