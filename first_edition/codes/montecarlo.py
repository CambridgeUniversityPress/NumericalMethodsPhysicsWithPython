# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from newtoncotes import f
import numpy as np

def montecarlo(f,a,b,n,option="uniform"):
    np.random.seed(314159)
    us = np.random.uniform(a, b, n)

    if option=="uniform":
        fs = f(us)
    else:
        c0 = 4 - 2*np.sqrt(2)
        c1 = -6 + 4*np.sqrt(2)
        xs = (-c0 + np.sqrt(2*c1*us + c0**2))/c1
        fs = f(xs)/(c0 + c1*xs)

    fbar, err = stats(fs)
    return (b-a)*fbar, (b-a)*err

def stats(fs):
    n = fs.size
    fbar = np.sum(fs)/n
    fsq = np.sum(fs**2)/n
    varfbar = (fsq - fbar**2)/(n - 1)
    return fbar, np.sqrt(varfbar)

if __name__ == '__main__':
    for n in 10**np.arange(2,7):
        avu, erru = montecarlo(f, 0., 1., n)
        avi, erri = montecarlo(f, 0., 1., n, option="is")
        rowf = "{0:7d}   {1:1.9f} {2:1.9f}   {3:1.9f} {4:1.9f}"
        print(rowf.format(n, avu, erru, avi, erri))
