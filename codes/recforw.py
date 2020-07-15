# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from math import exp

def forward(nmax=22):
    oldint = 1 - exp(-1)
    for n in range(1,nmax):
        print(n-1, oldint)
        newint = n*oldint - exp(-1)
        oldint = newint

print("n = 20 answer is 0.0183504676972562")
print("n, f[n]")
forward()
