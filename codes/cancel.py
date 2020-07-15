# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from math import exp, log

def f(x):
    return (exp(x) - 1)/x

def g(x):
    w = exp(x)
    if w==0.:
        val = -1/x
    elif w==1.:
        val = 1.
    else: 
        val = (w-1)/log(w)
    return val

xs = [10**(-i) for i in (14, 15, 16)] 
xs += [-10**(-i) for i in (15, 16, 17)]

fvals = [f(x) for x in xs] 
gvals = [g(x) for x in xs] 

for x, fval, gval in zip(xs, fvals, gvals):
    print(x, fval, gval)
