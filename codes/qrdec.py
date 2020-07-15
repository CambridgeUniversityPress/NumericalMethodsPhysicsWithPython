# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import testcreate
from power import mag
import numpy as np

def qrdec(A):
    n = A.shape[0]
    Ap = np.copy(A)
    Q = np.zeros((n,n))
    R = np.zeros((n,n))
    for j in range(n):
        for i in range(j):
            R[i,j] = Q[:,i]@A[:,j]
            Ap[:,j] -= R[i,j]*Q[:,i]

        R[j,j] = mag(Ap[:,j])
        Q[:,j] = Ap[:,j]/R[j,j]
    return Q, R

def testqrdec(A):
    n = A.shape[0]
    Q, R = qrdec(A)
    diffa = A - Q@R
    diffq = np.transpose(Q)@Q - np.identity(n) 
    print(n, mag(diffa), mag(diffq))

if __name__ == '__main__':
    for n in range(4,10,2):
        A, bs = testcreate(n,21)
        testqrdec(A)
