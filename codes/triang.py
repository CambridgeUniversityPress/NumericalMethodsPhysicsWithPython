# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def forsub(L,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in range(n):
        xs[i] = (bs[i] - L[i,:i]@xs[:i])/L[i,i]
    return xs

def backsub(U,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in reversed(range(n)):
        xs[i] = (bs[i] - U[i,i+1:]@xs[i+1:])/U[i,i]
    return xs

def testcreate(n,val):
    A = np.arange(val,val+n*n).reshape(n,n)
    A = np.sqrt(A)
    bs = (A[0,:])**2.1
    return A, bs

def testsolve(f,A,bs):
    xs = f(A,bs); print(xs)
    xs = np.linalg.solve(A,bs); print(xs)

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    L = np.tril(A)
    testsolve(forsub,L,bs)
    print(" ")
    U = np.triu(A)
    testsolve(backsub,U,bs)
