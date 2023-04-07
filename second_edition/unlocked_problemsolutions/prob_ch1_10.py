# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 1, problem 10

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np

# Note that the numpy functions are actually much more general,
# with optional arguments, functionality for higher-dimensional arrays, etc.
# Here we mock up the behavior for 1d arrays, so these functions no longer
# feel like black boxes

def myargmax(xs):
    val = xs[0]
    ival = 0
    for i,x in enumerate(xs):
        if x>val:
            val = x
            ival = i
    return ival

def mywhere(chs):
    for i,ch in enumerate(chs):
        if ch>0:
            val = np.array([i]),
            break
    else:
        val = np.array([]),
    return val

def myall(xs):
    val = True
    for x in xs:
        if x==False:
            val = False
    return val

xs = np.array([-0.2,-1.5,0.7,0.9,0.,0.5,-200.,1.1,-500.])
print(np.argmax(xs))
print(myargmax(xs))
print("")
print("")

xs = np.arange(10)*0.1
print(np.where(0.7==xs))
print(np.where(0.8==xs))
print(np.where(0.==xs))
print("")
print(mywhere(0.7==xs))
print(mywhere(0.8==xs))
print(mywhere(0.==xs))
print("")
print("")

xs = [True, False, False]
xs = [True, True, True]
xs = [1, 2, 0.]
xs = [1, 2, 3.]
print(np.all(xs))
print(myall(xs))

