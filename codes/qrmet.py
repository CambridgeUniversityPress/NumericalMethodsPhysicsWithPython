# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import testcreate
from qrdec import qrdec
import numpy as np

def qrmet(inA,kmax=100):
    A = np.copy(inA)
    for k in range(1,kmax):
        Q, R = qrdec(A)
        A = R@Q
        print(k, np.diag(A))

    qreigvals = np.diag(A)
    return qreigvals

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    qreigvals = qrmet(A,6)
    print(" ")
    npeigvals, npeigvecs = np.linalg.eig(A); print(npeigvals)
