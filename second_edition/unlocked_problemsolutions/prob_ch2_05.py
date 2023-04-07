# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 2, problem 5

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import log10
import matplotlib.pyplot as plt

def spike(x):
    return x**6 + log10(abs(1. + 3.*(1.-x)))/10.

def plot1():
    plt.xlabel('x')
    plt.ylabel('spike(x)')

    xs = [0.01*i + 0.5 for i in range(100)]
    fs = [spike(x) for x in xs]
    plt.plot(xs,fs,'r-')

    plt.show()

def plot2():
    plt.xlabel('x')
    plt.ylabel('spike(x)')

    xs = [0.000005*i + 1.33 for i in range(1000)]
    fs = [spike(x) for x in xs]
    plt.plot(xs,fs,'r-')

    plt.show()


plot1()
plot2()

# The lesson here is that if you use too coarse a grid (and don't think about
# the analytic properties of your function) you might simply miss the spike.
