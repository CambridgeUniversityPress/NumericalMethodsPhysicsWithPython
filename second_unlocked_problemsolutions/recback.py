# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from math import exp

def backward(nmax=31):
    oldint = 0.01
    for n in reversed(range(20,nmax)):
        print(n, oldint)
        newint = (oldint + exp(-1))/n
        oldint = newint

print("n = 20 answer is 0.0183504676972562")
print("n, f[n]")
backward()
