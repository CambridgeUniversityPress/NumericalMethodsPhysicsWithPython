# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from newtoncotes import f, rectangle, trapezoid, simpson
import numpy as np

def adaptive(f,a,b,integrator,kmax = 20,tol = 1.e-12):
    functodenom = {rectangle:1, trapezoid:3, simpson:15}
    denom = functodenom[integrator]

    n = 2
    val = integrator(f,a,b,n) 
    for k in range(kmax):
        nprime = 2*n-1
        valprime = integrator(f,a,b,nprime)
        err = abs(valprime-val)/denom
        err /= abs(valprime)
        print(nprime, valprime, err)
        if err<tol:
            break

        n, val = nprime, valprime
    else:
        valprime = None
    return valprime

if __name__ == '__main__':
    ans = np.log(1 + np.sqrt(2))
    print(ans); print("")
    for integrator in (rectangle, trapezoid, simpson):
        print(adaptive(f, 0., 1., integrator)); print("")
