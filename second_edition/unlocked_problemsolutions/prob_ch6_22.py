# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 6, problem 22

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import exp, sin, cos, log10, pi
import matplotlib.pyplot as plt
import numpy as np
from triginterp import generatedata, computeparams, triginterp

# this is the new function to investigate.
def f(x):
    return np.exp(np.sin(x) + np.cos(x))

# this merely does the plotting. Note that all the heavy lifting
# is done in the functions we imported from triginterp.py
def plotall(xs,ys,interpxs,interpys):
    plt.xlabel('$x$', fontsize=18)
    plt.plot(interpxs,f(interpxs),'b-',label='$f(x)$',linewidth=1.5)
    plt.plot(xs,ys,'ro',label='$(x_j, y_j)$')
    plt.plot(interpxs,interpys,'g--',label='$p(x)$')
    plt.xlim(-0.25,2*pi+0.25)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    # be sure to distinguish between x,y for the data and x,y for plotting purposes
    # note that this will generate two separate plots (which is what the question is asking)
    for n in 6, 8:
        dataxs, datays = generatedata(n, f)
        aparams, bparams = computeparams(dataxs, datays)
        x = 0.3; pofx = triginterp(aparams, bparams, x)
        interpxs = (2*pi/500)*np.arange(501)
        interpys = triginterp(aparams,bparams,interpxs)
        plotall(dataxs,datays,interpxs,interpys)

