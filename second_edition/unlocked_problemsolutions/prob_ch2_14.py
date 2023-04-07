# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 2, problem 14

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import pi, sin, factorial

# n-th term is (-1)**n * x**(2*n+1)/(2*n+1)!
# factor relating n-th and (n-1)-th term is -x**2/((2*n+1)*2*n)

def compsin(x):
    n = 0
    oldsum, newsum, term = 0., x, x
    while newsum != oldsum:
        oldsum = newsum
        n += 1
        term *= -x**2/((2*n+1)*2*n)
        newsum += term
        #print(n, newsum, term)
    return newsum

x=0.1
print(x, 'library sin', sin(x))
print(x, 'compsin', compsin(x))

x=40.
print(x, 'library sin', sin(x))
print(x, 'compsin', compsin(x))
print(x, 'compsin smarter', compsin(x - 12*pi)) # uses sin(theta + k*2pi) = sin(theta)
# this uses a trigonometric identity to go to a smaller value of the argument
