# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 1, problem 3

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

xs = [1.3, 4.5, -2.4, -7.5]

# Pythonic, but too long
for i, x in reversed(list(enumerate(xs))):
    print(i, x)

print("")

# unPythonic
n = len(xs)
for i in reversed(range(n)):
    print(i, xs[i])
