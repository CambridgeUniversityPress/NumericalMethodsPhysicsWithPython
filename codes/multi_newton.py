# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from gauelim_pivot import gauelim_pivot
from jacobi import termcrit
import numpy as np

def fs(xs):
    x0, x1 = xs
    f0 = x0**2 - 2*x0 + x1**4 - 2*x1**2 + x1
    f1 = x0**2 + x0 + 2*x1**3 - 2*x1**2 - 1.5*x1 - 0.05
    return np.array([f0,f1])

def jacobian(fs,xs,h=1.e-4):
    n = xs.size
    iden = np.identity(n)
    Jf = np.zeros((n,n))
    fs0 = fs(xs)
    for j in range(n):
        fs1 = fs(xs+iden[:,j]*h)
        Jf[:,j] = (fs1 - fs0)/h
    return Jf, fs0

def multi_newton(fs,jacobian,xolds,kmax=200,tol=1.e-8):
    for k in range(1,kmax):
        Jf, fs_xolds = jacobian(fs, xolds)
        xnews = xolds + gauelim_pivot(Jf, -fs_xolds)

        err = termcrit(xolds,xnews)
        print(k, xnews, err)
        if err < tol:
            break

        xolds = np.copy(xnews)
    else:
        xnews = None
    return xnews

if __name__ == '__main__':
    xolds = np.array([1.,1.])
    xnews = multi_newton(fs,jacobian,xolds)
    print(xnews); print(fs(xnews))
