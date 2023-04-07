# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 5, problem 11

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import exp, sqrt, log

# the following two functions correspond to the initial problem
# before root suppression. You don't actually need them for this problem
def f(x):
    return exp(x - sqrt(x)) - x

def fprime(x):
    return -1 + exp(x - sqrt(x))*(1 - 1/(2*sqrt(x)))

# the next two functions correspond to the root-suppression problem
def u(x):
    return (exp(x - sqrt(x)) - x)/(x-1)

# Newton's method needs the derivative as input
# here we take the derivative analytically (even though it's somewhat messy)
def uprime(x):
    return -(exp(x - sqrt(x)) - x)/(x-1)**2 +  (-1 + exp(x - sqrt(x))*(1 - 1/(2*sqrt(x))))/(x-1)

# the previous two functions are specific to this problem, but the following one is general
def newton(f,fprime,xold,kmax=200,eps=1.e-8):
    for k in range(1,kmax):
        xnew = xold - f(xold)/fprime(xold)
        xdiff = xnew - xold
        print("{0:2d} {1:1.16f} {2:1.16f}".format(k,xnew,xdiff))

        if abs(xdiff/xnew) < eps:
            break

        xold = xnew
    else:
        xnew = None

    return xnew

if __name__ == '__main__':
    # no matter where you start, you always suppress the root you already know
    root = newton(u,uprime,2.0); print(root)
    root = newton(u,uprime,0.5); print(root)
    root = newton(u,uprime,4.0); print(root)
    root = newton(u,uprime,0.1); print(root)

    print("")
    # you don't actually need to do the following comparison, but you may 
    # be interested in seeing how Newton's method behaves without root suppression.
    root = newton(f,fprime,2.0); print(root)
    root = newton(f,fprime,0.5); print(root)
    root = newton(f,fprime,4.0); print(root)
    root = newton(f,fprime,0.1); print(root)
