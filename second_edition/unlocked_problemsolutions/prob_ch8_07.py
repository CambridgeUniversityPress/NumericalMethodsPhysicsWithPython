# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 8, problem 7

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ivp_one import f, euler

# SciPy's odeint() uses the opposite order for the arguments
# in comparison to our own convention.
def g(y,x):
    return ((2*x)/(1-x**2))*y - (30/(1-x**2)) - y**2

def adams(f,a,b,n,yinit):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros(n)

    y0 = yinit
    ys[0] = y0
    # This is the point where we approximate the second starting value
    y1 = y0 + h*f(xs[0],y0)
    for j,x in enumerate(xs[:-1]):
        ys[j+1] = y1
        # This is a straightforward implementation of the Adams-Bashforth prescription
        y2 = y1 + (h/2)*(3*f(xs[j+1],y1) - f(xs[j],y0))
        y0, y1 = y1, y2
    return xs, ys

# We are comparing against Euler (which we employ to get started)
# as well as the "exact" solution from SciPy.
def plot_full():
    plt.xlabel('$x$', fontsize=18)
    plt.ylabel('$y(x)$', fontsize=18)

    a,b,n,y0 = 0.05, 0.49, 22, 19.53
    xs, ys = adams(f,a,b,n,y0)  
    plt.plot(xs, ys, 'bs', label='Adams', linewidth=2)

    xs, ys = euler(f,a,b,n,y0)  
    print(xs)
    print(ys)
    plt.plot(xs, ys, 'r^', label='Euler', linewidth=2)

    xs = np.linspace(a,b,1000)
    ys = odeint(g,y0,xs)  
    plt.plot(xs, ys, 'y-', label='Exact', linewidth=2)

    plt.legend(loc=1)

    plt.xlim(a-0.05, b+0.05)
    plt.ylim(-55,25)

    plt.show()
        
if __name__ == '__main__':
    plot_full()
    # The results are fascinating: even though our second point (coming
    # from Euler) is bad, from the third point onwards the Adams-Bashforth
    # answer is considerably closer to the true answer than the Euler one was.
