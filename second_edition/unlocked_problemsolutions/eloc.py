# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

import numpy as np

def psi(al, oms, X):
    rexp = np.sum(oms*X**2)
    return np.exp(-0.5*al*rexp)

def ekin(al, oms, X, h=0.01, hom = 1.):
    npart, ndim = X.shape
    psiold = psi(al, oms, X)
    kin = 0.
    for (j,el), r in np.ndenumerate(X):
        X[j,el] = r + h
        psip = psi(al, oms, X)
        X[j,el] = r - h
        psim = psi(al, oms, X)
        X[j,el] = r
        lapl = (psip + psim - 2.*psiold)/h**2
        kin += -0.5*hom*lapl/psiold
    return kin

def epot(oms, X, strength = 3, m = 1., q = 1.):
    npart, ndim = X.shape
    pot = 0.5*m*np.sum(oms**2*X**2)
    for k in range(1,npart):
        for j in range(k):
            r2 = np.sum((X[j,:] - X[k,:])**2)
            pot += strength*np.exp(-q*r2)
    return pot

if __name__ == '__main__':
    npart, ndim, al = 4, 3, 0.6
    oms = np.arange(1, 1 + ndim)
    X = 1/np.arange(1, 1 + npart*ndim).reshape(npart, ndim)
    print(ekin(al, oms, X), epot(oms, X))
