# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from math import exp, sqrt

def g(x):
    return exp(x - sqrt(x))

def fixedpoint(g,xold,kmax=200,tol=1.e-8):
    for k in range(1,kmax):
        xnew = g(xold)

        xdiff = xnew - xold
        print("{0:2d} {1:1.16f} {2:1.16f}".format(k,xnew,xdiff))

        if abs(xdiff/xnew) < tol:
            break

        xold = xnew
    else:
        xnew = None

    return xnew

if __name__ == '__main__':
    for xold in (0.99, 2.499):
        x = fixedpoint(g,xold)
        print(x)
