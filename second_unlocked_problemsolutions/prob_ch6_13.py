# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 6, problem 13

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi
from barycentric import f, generatedata

# indeed, this is a simplified version of splineinterp() from splines.py
def linsplin(dataxs,datays,x):
    k = np.argmax(dataxs>x)
    xk = dataxs[k]; xk1 = dataxs[k-1]
    yk = datays[k]; yk1 = datays[k-1]

    val = yk1*(xk-x)/(xk-xk1) + yk*(x-xk1)/(xk-xk1)
    return val

# this is a quick-and-dirty implementation. One could, instead,
# pass in n as an argument and then call the present function twice.
def plotting():
    plt.xlabel('$x$', fontsize=24)

    n = 15
    dataxs, datays = generatedata(n, f, "equi")
    interpxsA = np.linspace(-1,1,100*n)
    interpysA = [linsplin(dataxs, datays, x) for x in interpxsA] 
    plt.plot(interpxsA, interpysA, 'b--', label='Splines $n=15$', linewidth=3)

    n = 150
    dataxs, datays = generatedata(n, f, "equi")
    interpxsB = np.linspace(-1,1,100*n)
    interpysB = [linsplin(dataxs, datays, x) for x in interpxsB] 
    plt.plot(interpxsB, interpysB, 'g-.', label='Splines $n=150$', linewidth=4)

    plt.plot(interpxsB, f(interpxsB), 'k-', label='Runge function', linewidth=2.5)

    plt.legend(loc=1)

    plt.xlim(-1.05, 1.05)
    plt.ylim(-0.15, 1.45)
    plt.show()

if __name__ == '__main__':
    plotting()
