# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 7, problem 21

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

# while it looks like we have to import too many things
from legendre import legendre
from legroots import legroots
from gauleg import gauleg_params
from newtoncotes import f
import numpy as np
import matplotlib.pyplot as plt

# this pays off in the main program, where we simply use
# the functionality we've already developed, basically just plotting.
if __name__ == '__main__':
    for n in 5,50,100:
        plt.xlabel('$position$', fontsize=25)
        plt.ylabel('$weight$', fontsize=25)
        xs, cs = gauleg_params(n)
        plt.plot(xs, cs, 'ro')
        plt.show()

# As mentioned in section 5.3.2, the x's obey a clustering property 
# around the ends of the interval. Note the scale on the y-axis:
# as you make n larger, the weights get smaller.
