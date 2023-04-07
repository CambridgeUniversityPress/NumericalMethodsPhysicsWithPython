# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 4, problem 24

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from triang import testcreate
import numpy as np

def mag(xs):
    return np.sqrt(np.sum(xs*xs))

# this is a modified version of the function from power.py
def power(A,kmax=6):
    zs = np.ones(A.shape[0])
    qs = zs/mag(zs)
    pn = p1 = 0. # this is the assumption to get things going
    for k in range(1,kmax):
        zs = A@qs
        qs = zs/mag(zs)
        p2 = qs@A@qs
        ait = pn - (p1 - pn)**2/(p2 - 2.*p1 + pn) # this is the core Aitken extrapolation 
        print(k,p2,ait)
        pn = p1
        p1 = p2

    lam = qs@A@qs
    return lam, qs

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    print(np.linalg.eigvals(A)[0])
    power(A)
