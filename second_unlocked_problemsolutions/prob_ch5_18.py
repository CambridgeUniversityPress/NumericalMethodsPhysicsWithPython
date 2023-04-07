# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 5, problem 18

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import cos, pi, sqrt
import numpy as np

# this is taken from psis.py
def hermite(n,x):
    if n==0:
        val2 = 1.
        dval2 = 0.
    elif n==1:
        val2 = 2*x
        dval2 = 2.
    else:
        val0 = 1.; val1 = 2*x
        for j in range(1,n):
            val2 = 2*x*val1 - 2*j*val0
            val0, val1 = val1, val2
        dval2 = 2*n*val0
    return val2, dval2

# this is very similar to legroots.py
def hermnewton(n,xold,kmax=200,tol=1.e-8):
    for k in range(1,kmax):
        val, dval = hermite(n,xold)
        xnew = xold - val/dval

        xdiff = xnew - xold
        if abs(xdiff/xnew) < tol:
            break

        xold = xnew
    else:
        xnew = None
    return xnew

# The only subtlety is in the line giving nume,
# which you should also try to improve upon.
# The rest looks a lot like legroots.py
def hermroots(n):
    roots = [0.]*n
    npos = n//2
    for i in range(npos):
        nume = 3.5*i+2 if n%2==0 else 3.6*i+3.6
        xold = nume/sqrt(2*n+1)
        root = hermnewton(n,xold) 
        roots[i] = -root
        roots[-1-i] = root
    return roots

# For not-too-large values of n, we agree with the library output
if __name__ == '__main__':
    for n in range(2,21):
        xs, ws = np.polynomial.hermite.hermgauss(n)
        roots = hermroots(n)
        print(n)
        print(np.abs(np.sort(roots)-np.sort(xs))<0.01)
        print("")
