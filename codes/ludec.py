# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import forsub, backsub, testcreate, testsolve
import numpy as np

def ludec(A):
    n = A.shape[0]
    U = np.copy(A)
    L = np.identity(n)

    for j in range(n-1):
        for i in range(j+1,n):
            coeff = U[i,j]/U[j,j]
            U[i,j:] -= coeff*U[j,j:]
            L[i,j] = coeff

    return L, U

def lusolve(A,bs):
    L, U = ludec(A)
    ys = forsub(L,bs)
    xs = backsub(U,ys)
    return xs

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    testsolve(lusolve,A,bs)
