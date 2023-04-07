# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from math import exp

def compexp(x):
    n = 0
    oldsum, newsum, term = 0., 1., 1.
    while newsum != oldsum:
        oldsum = newsum
        n += 1
        term *= x/n
        newsum += term
        print(n, newsum, term)
    return newsum

for x in (0.1, 20., -20.):
    print("x, library exp(x):", x, exp(x))
    val = compexp(x)
