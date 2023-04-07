# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 5, problem 21

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from gauelim_pivot import gauelim_pivot
from jacobi import termcrit
import numpy as np
from math import exp, sqrt

# Only subtlety is that this is g (not f)
# that is, we have moved x0 (in the first equation) and x1 (in the 
# second equation) to the left-hand side.
def gs(xs):
    x0, x1 = xs
    g0 = (x0**2 - 2*x1**2 + x1)/2
    g1 = (x0**2 + x0 - 2*x1**2 - 0.05)/1.5
    return np.array([g0,g1])

# This is a general implementation.
def fixedsystem(gs,xolds,kmax=100,tol=1.e-8):
    for k in range(1,kmax):
        xnews = gs(xolds)

        xdiffs = xnews - xolds
        print(k,xnews,xdiffs)

        if abs(np.sum(xdiffs/xnews)) < tol:
            break

        xolds = np.copy(xnews)
    else:
        xnews = None

    return xnews

if __name__ == '__main__':
    # Just like in the one-dimensional case, you keep
    # finding one of the solutions (or you crash).
    # Again, like in 1d, the g we picked is not unique, 
    # so you may have better luck with the other solution
    # if you rewrite the starting problem differently.
    # You should try out other initial guesses as well.
    xolds = np.array([1., 1.])
    #xolds = np.array([-0.35, -0.45])
    #xolds = np.array([1., -1.])
    xnews = fixedsystem(gs,xolds)
    print(xnews)
