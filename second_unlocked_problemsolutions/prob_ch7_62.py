# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 7, problem 62

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np
from montecarlo import montecarlo

def f(x, A=0.25):
    return np.exp(-x**2)*np.exp(A*x)

def g(x, A):
    return np.exp(A*x)

# The derivation is straightforward, but you should use
# paper and pencil to do this. Note that the shift 
# the question asks you to do relies on employing
# an infinite interval. 
# We evaluate this here as a benchmark for what follows:
A = 0.25
print("analytical answer: ", np.sqrt(np.pi)*np.exp(A**2/4))

# Make sure you use enough points so the sample error is small
nsteps = 10**8
avu, erru = montecarlo(f,-10,10, nsteps)
rowf = "{0:7d}   {1:1.9f} {2:1.9f}"
print(rowf.format(nsteps, avu, erru))

# you should try to do the same thing using Problem 7.25
xs, ws = np.polynomial.hermite.hermgauss(2)
print("n=2: ", np.sum(ws*g(xs,A)))
xs, ws = np.polynomial.hermite.hermgauss(3)
print("n=3: ", np.sum(ws*g(xs,A)))
# It is striking how well we can do with merely 3 points

