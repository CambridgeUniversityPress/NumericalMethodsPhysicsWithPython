# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 5, problem 15

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import exp, sqrt, log
import matplotlib.pyplot as plt

# our example function
def f(x):
    return exp(x - sqrt(x)) - x

# a general implementation of the Steffensen method
# note that no derivative input is needed
def steffensen(f,xold,kmax=200,eps=1.e-8):
    for k in range(1,kmax):
        g = (f(xold + f(xold)) - f(xold))/f(xold)
        xnew = xold - f(xold)/g
        xdiff = xnew - xold
        print("{0:2d} {1:1.16f} {2:1.16f}".format(k,xnew,xdiff))

        if abs(xdiff/xnew) < eps:
            break

        xold = xnew
    else:
        xnew = None

    return xnew

if __name__ == '__main__':
    # depending on the starting guess, we find one or the other root
    root = steffensen(f,3.0); print(root)
    root = steffensen(f,0.5); print(root)
