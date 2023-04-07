# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 2, problem 19

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import scipy.special as sp
import matplotlib.pyplot as plt

# Note that this function (from the main text) is wasteful, 
# since it evaluates all the polynomials up to n also, but throws them away
def legendre(n,x):
    if n==0:
        val2 = 1.
        dval2 = 2.
    elif n==1:
        val2 = x
        dval2 = 1.
    else:
        val0 = 1.; val1 = x
        for i in range(1,n):
            val2 = ((2*i+1)*x*val1 - i*val0)/(i+1)
            val0, val1 = val1, val2
        dval2 = n*(val0-x*val1)/(1.-x**2)
    return val2, dval2

# This, too, is wasteful.
def legendre_other(n,x):
    if n==0:
        val2 = 1.
    elif n==1:
        val2 = x
    else:
        val0 = 1.; val1 = x
        for i in range(1,n):
            val2 = 2*x*val1 - val0  - (x*val1 - val0)/(i+1)
            val0, val1 = val1, val2
    return val2

def testlegendre(n):
    xs = [0.01*i for i in range(100)]
    for x in xs:
        print(x,legendre(n,x)[0], legendre_other(n,x), sp.eval_legendre(n,x))

def plotlegendre(n, der=0, nsteps=1000):
    plt.xlabel('$x$', fontsize=20)

    dertostr = {0: "$P_n(x)$", 1: "$P_n'(x)$"}
    plt.ylabel("$P^A_n(x) - P^A_n(x)$", fontsize=20)
        
    xs = [i/nsteps for i in range (-nsteps+1,nsteps)]
    ys = [legendre(n,x)[der]-legendre_other(n,x) for x in xs]
    plt.plot(xs, ys, 'ro', linewidth=3)

    plt.show()

if __name__ == '__main__':
    # The question doesn't specify how to do the comparison,
    # so we provide both a table and a plot
    testlegendre(1000)
    plotlegendre(1000)
