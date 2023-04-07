# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 6, problem 5

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np
import matplotlib.pyplot as plt
from barycentric import weights, bary

# We use different names/functions for the two parts of the problem
def f(x):
    return np.exp(np.sin(20*x))

def g(x):
    return 100*(x-1)**3*np.cos(4*(x-1))*np.exp(5*(x-1))

# Only subtlety is that we use one plotting function to solve both parts,
# so we pass in as arguments the name of the function and the number of points
def plotting(func, ns):
    plt.xlabel('$x$', fontsize=18)

    n = ns[0]
    dataxs = -np.cos(np.linspace(0,np.pi,n))
    datays = func(dataxs)
    interpxsA = np.linspace(-1,1,100*n)
    ws = weights(dataxs)
    interpysA = [bary(dataxs, datays, ws, x) for x in interpxsA] 
    plt.plot(interpxsA, interpysA, 'b--', label='Lagrange $n={0}$'.format(n), linewidth=2)

    n = ns[1]
    dataxs = -np.cos(np.linspace(0,np.pi,n))
    datays = func(dataxs)
    interpxsB = np.linspace(-1,1,100*n)
    ws = weights(dataxs)
    interpysB = [bary(dataxs, datays, ws, x) for x in interpxsB] 
    plt.plot(interpxsB, func(interpxsB), 'k-', label='Underlying function', linewidth=1.5)
    plt.plot(interpxsB, interpysB, 'g-.', label='Lagrange $n={0}$'.format(n), linewidth=3)

    plt.legend(loc=1)
    plt.xlim(-1.05, 1.05)
    plt.show()

if __name__ == '__main__':
    plotting(f, (15,60))
    plotting(g, (7,15))
