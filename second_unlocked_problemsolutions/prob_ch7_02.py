# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 7, problem 2

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from newtoncotes import f, rectangle, trapezoid, simpson
import numpy as np

# this is very similar to the newtoncotes functions
# the main/only difference being the evaluation at xs + h/2
def midpoint(f,a,b,n):
    h = (b-a)/(n-1)
    xs = a + np.arange(n-1)*h
    contribs = f(xs + h/2)
    return h*np.sum(contribs)

# The analytic derivation is very similar to Eq. (7.48) and (7.49), but you
# have to be careful with factors of 2. It leads to Eq. (7.21). Then, by
# an argument like that in Eq. (7.18) you get Eq. (7.22).

if __name__ == '__main__':
    ans = np.log(1 + np.sqrt(2))
    print(ans)

    for integrator in (rectangle, trapezoid, midpoint, simpson):
        print(integrator(f, 0., 1., 51), end=" ")
    print("")

# We get 5 decimal digits of agreements, just like in the trapezoid rule
# (but much better than rectangle rule). Obviously, still worse than simpson's rule.

