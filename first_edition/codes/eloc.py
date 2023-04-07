# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def psi(al, oms, rs):
    rexp = np.sum(oms*rs**2)
    return np.exp(-0.5*al*rexp)

def ekin(al, oms, rs, h=0.01, hom = 1.):
    npart, ndim = rs.shape
    psiold = psi(al, oms, rs)
    kin = 0.
    for j in range(npart):
        numer = 0.
        for el in range(ndim):
            r = rs[j,el]
            rs[j,el] = r + h
            psip = psi(al, oms, rs)
            rs[j,el] = r - h
            psim = psi(al, oms, rs)
            rs[j,el] = r
            numer += psip + psim - 2.*psiold
        lapl = numer/h**2
        kin += -0.5*hom*lapl/psiold
    return kin

def epot(oms, rs, strength = 3, m = 1.):
    npart, ndim = rs.shape
    pot = 0.5*m*np.sum(oms**2*rs**2)
    for k in range(1,npart):
        for j in range(k):
            r2 = np.sum((rs[j,:] - rs[k,:])**2)
            pot += strength*np.exp(-r2)
    return pot

if __name__ == '__main__':
    npart, ndim, al = 4, 3, 0.6
    oms = np.arange(1, 1 + ndim)
    rs = np.arange(1, 1 + npart*ndim).reshape(npart, ndim)
    rs = 1/rs
    print(ekin(al, oms, rs), epot(oms, rs))
