# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 8, problem 29

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from ivp_two import rk4_gen
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Right-hand side of our ODE for Numerov
def f(x):
    return 5/(1+x**2)

# Right-hand side of our ODE for the Euler/RK2/RK4/etc approaches
def fsforrk4(x,yvals):
    y0, y1 = yvals
    f0 = y1
    f1 = 5/(1+x**2)*y0
    return np.array([f0, f1])

# This is similar to our other IVP solvers. The main difference
# is that ys here is a 1d array (even though we are tackling a 2nd-order ODE).    
def numerov(f,a,b,n,yinits):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros(n)

    y0 = yinits[0]
    ys[0] = y0
    # This is the point where we approximate the second starting value
    y1 = y0 + h*yinits[1]
    for j,x in enumerate(xs[:-2]):
        ys[j+1] = y1
        # This is a straightforward implementation of the Numerov prescription
        y2 = (2 + 5*h**2*f(xs[j+1])/6)*y1 - (1 - h**2*f(xs[j])/12)*y0
        y2 /= 1 - h**2*f(xs[j+2])/12
        y0, y1 = y1, y2
    ys[-1] = y1
    return xs, ys
        
# The problem doesn't ask us to plot the two solutions, but we do so anyway.
if __name__ == '__main__':
    a, b, n = 0., 1, 500
    yinits = np.array([1.,0])
    xs, ys = numerov(f,a,b,n,yinits)
    plt.plot(xs, ys, 'rD', label='Numerov', linewidth=2)

    yinits = np.array([1.,0])
    xs, ys = rk4_gen(fsforrk4,a,b,n,yinits)
    plt.plot(xs, ys[:,0], 'bo', label='RK4', linewidth=2)

    plt.legend(loc=2)

    plt.show()
