# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import backsub, testcreate, testsolve
import numpy as np

def gauelim_pivot(inA,inbs):
    A = np.copy(inA)
    bs = np.copy(inbs)
    n = bs.size

    for j in range(n-1):
        k = np.argmax(np.abs(A[j:,j])) + j
        if k != j:
            A[j,:], A[k,:] = A[k,:], A[j,:].copy()
            bs[j], bs[k] = bs[k], bs[j]

        for i in range(j+1,n):
            coeff = A[i,j]/A[j,j]
            A[i,j:] -= coeff*A[j,j:]
            bs[i] -= coeff*bs[j]

    xs = backsub(A,bs)
    return xs

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    testsolve(gauelim_pivot,A,bs)
