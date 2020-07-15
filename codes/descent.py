# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from jacobi import termcrit
import numpy as np

def phi(xs):
    x0, x1 = xs
    return x0**2 - 2*x0 + x1**4 - 2*x1**2 + x1

def gradient(phi,xs,h=1.e-6):
    n = xs.size
    phi0 = phi(xs)
    Xph = (xs*np.ones((n,n))).T + np.identity(n)*h
    grad = (phi(Xph) - phi0)/h
    return grad

def descent(phi,gradient,xolds,gamma=0.15,kmax=200,tol=1.e-8):
    for k in range(1,kmax):
        xnews = xolds - gamma*gradient(phi,xolds)

        err = termcrit(xolds,xnews)
        print(k, xnews, err, phi(xnews))
        if err < tol:
            break

        xolds = np.copy(xnews)
    else:
        xnews = None
    return xnews

if __name__ == '__main__':
    xolds = np.array([2.,0.25])
    xnews = descent(phi, gradient, xolds)
    print(xnews)
