# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from bisection import f as phi
from math import sqrt

def golden(phi, x0, x1, kmax=200, tol=1.e-8):
    varphi = 0.5*(1 + sqrt(5))
    for k in range(1,kmax):
        x2 = x0 + (x1-x0)/(varphi+1)
        x3 = x0 + (x1-x0)*varphi/(varphi+1)
        
        if phi(x3) < phi(x2):
            x0 = x2
        else:
            x1 = x3
            
        xnew = (x0+x1)/2
        xdiff = abs(x1-x0)
        rowf = "{0:2d} {1:1.16f} {2:1.16f} {3:1.16f}"
        print(rowf.format(k, xnew, xdiff, phi(xnew)))

        if abs(xdiff) < tol:
            break
    else:
        xnew = None
    return xnew

if __name__ == '__main__':
    val = golden(phi,0.,3.5); print(val)
