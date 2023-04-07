# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 8, problem 28

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from ivp_two import fs, rk4_gen
import numpy as np
import matplotlib.pyplot as plt

# The interface is identical to many of our other functions. A subtlety is noted below.
def verlet(f,a,b,n,yinit):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros((n, yinits.size))

    yvals = np.copy(yinits)
    for j,x in enumerate(xs):
        ys[j,:] = yvals
        k = yvals[1] + (h/2)*f(x,yvals)[1]
        # Note that the following two updates need to be done in the right order
        yvals[0] += h*k
        yvals[1] = k + (h/2)*f(x+h,[yvals[0],k])[1] # this uses the previous line
    return xs, ys

if __name__ == '__main__':
    # n=6 is artificially small, precisely in order to bring out the differences
    # in the results coming from the two methods.
    for n in 6, 100:
        a, b = 0.05, 0.49
        yinits = np.array([0.0926587109375, 1.80962109375])

        xs, ys = rk4_gen(fs,a,b,n,yinits)
        plt.plot(xs, ys[:,0], 'ro', label='RK4', linewidth=2)

        xs, ys = verlet(fs,a,b,n,yinits)
        plt.plot(xs, ys[:,0], 'bD', label='Verlet', linewidth=2)

        plt.legend()
        plt.show()
