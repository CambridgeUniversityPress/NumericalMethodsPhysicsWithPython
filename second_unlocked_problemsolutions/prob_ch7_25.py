# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 7, problem 25

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from prob_ch5_18 import hermroots # solution to Problem 5.18
from psis import hermite
import numpy as np
from math import factorial

# Straigthforward implementations, similar to gauleg.py
def gauherm_params(n):
    xs = np.array(hermroots(n))
    cs = 2**(n+1)*np.sqrt(np.pi)*factorial(n)/(hermite(n,xs)[1]**2)
    return xs, cs

def gauherm(f,n):
    xs, cs = gauherm_params(n)
    return np.sum(cs*f(xs))

if __name__ == '__main__':
    n = 20
    xs, cs = np.polynomial.hermite.hermgauss(n)
    print(xs)
    print(cs)
    print("")
    xs, cs = gauherm_params(n)

    # You should spend some time thinking about why the
    # the sorting is necessary here
    print(np.sort(xs))
    # Similarly, you should also think why argsort is needed here.
    inds = np.argsort(xs)
    print(cs[inds])
