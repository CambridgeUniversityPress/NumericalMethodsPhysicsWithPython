# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import testcreate
from invpowershift import invpowershift
from qrmet import qrmet
import numpy as np

def eig(A,eps=1.e-12):
    n = A.shape[0]
    eigvals = np.zeros(n)
    eigvecs = np.zeros((n,n))
    qreigvals = qrmet(A)
    for i, qre in enumerate(qreigvals):
        eigvals[i], eigvecs[:,i] = invpowershift(A,qre+eps)
    return eigvals, eigvecs

def testeigall(f,A):
    eigvals, eigvecs = f(A)
    npeigvals, npeigvecs = np.linalg.eig(A)
    print(eigvals); print(npeigvals)
    print(" ")
    for eigvec, npeigvec in zip(eigvecs.T,npeigvecs.T):
        print(eigvec); print(npeigvec)
        print(" ")
    
if __name__ == '__main__':
    A, bs = testcreate(4,21)
    testeigall(eig,A)
