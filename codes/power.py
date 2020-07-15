# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import testcreate
import numpy as np

def mag(xs):
    return np.sqrt(np.sum(xs*xs))

def power(A,kmax=6):
    zs = np.ones(A.shape[0])
    qs = zs/mag(zs)
    for k in range(1,kmax):
        zs = A@qs
        qs = zs/mag(zs)
        print(k,qs)

    lam = qs@A@qs
    return lam, qs

def testeigone(f,A,indx=0):
    eigval, eigvec = f(A)
    print(" "); print(eigval); print(eigvec)
    npeigvals, npeigvecs = np.linalg.eig(A)
    print(" ")
    print(npeigvals[indx]); print(npeigvecs[:,indx])

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    testeigone(power,A)
