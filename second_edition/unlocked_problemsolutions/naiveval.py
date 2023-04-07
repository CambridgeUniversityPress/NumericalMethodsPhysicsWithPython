# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from math import sqrt

def naiveval(x):
    return 1/(sqrt(x**2 + 1) - x)

xs = [10**i for i in range(4,8)] 
ys = [naiveval(x) for x in xs] 
for x, y in zip(xs, ys):
    print(x, y)
