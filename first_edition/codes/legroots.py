# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from legendre import legendre
import numpy as np

def legnewton(n,xold,kmax=200,tol=1.e-8):
    for k in range(1,kmax):
        val, dval = legendre(n,xold)
        xnew = xold - val/dval

        xdiff = xnew - xold
        if abs(xdiff/xnew) < tol:
            break

        xold = xnew
    else:
        xnew = None
    return xnew

def legroots(n):
    roots = np.zeros(n)
    npos = n//2
    for i in range(npos):
        xold = np.cos(np.pi*(4*i+3)/(4*n+2))
        root = legnewton(n,xold) 
        roots[i] = -root
        roots[-1-i] = root
    return roots

if __name__ == '__main__':
    roots = legroots(9); print(roots)
