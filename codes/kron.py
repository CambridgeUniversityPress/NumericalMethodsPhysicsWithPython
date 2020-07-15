# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def paulimatrices():
    sigx = np.array([0.,1,1,0]).reshape(2,2)
    sigy = np.array([0.,-1j,1j,0]).reshape(2,2)
    sigz = np.array([1.,0,0,-1]).reshape(2,2)
    return sigx, sigy, sigz

def kron(U,V):
    n = U.shape[0]
    p = V.shape[0]
    W = np.zeros((n*p,n*p), dtype=np.complex64)
    for i in range(n):
        for k in range(n):
            for j in range(p):
                for l in range(p):
                    W[p*i+j,p*k+l] = U[i,k]*V[j,l]
    return W

if __name__ == '__main__':
    sigx, sigy, sigz = paulimatrices()
    allones = np.ones((3,3))
    kronprod = kron(sigx,allones); print(kronprod.real)
