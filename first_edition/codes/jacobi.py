# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import testcreate, testsolve
import numpy as np

def termcrit(xolds,xnews):
    errs = np.abs((xnews - xolds)/xnews)
    return np.sum(errs)

def jacobi(A,bs,kmax=50,tol=1.e-6):
    n = bs.size
    xnews = np.zeros(n)

    for k in range(1,kmax):
        xs = np.copy(xnews)

        for i in range(n):
            slt = A[i,:i]@xs[:i]
            sgt = A[i,i+1:]@xs[i+1:]
            xnews[i] = (bs[i] - slt - sgt)/A[i,i]

        err = termcrit(xs,xnews)
        print(k, xnews, err)
        if err < tol:
            break
    else:
        xnews = None

    return xnews

if __name__ == '__main__':
    n = 4; val = 21
    A, bs = testcreate(n,val)
    A += val*np.identity(n)
    testsolve(jacobi,A,bs)
