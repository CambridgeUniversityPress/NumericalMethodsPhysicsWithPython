# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 4, problem 6

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np

# The next two functions are the same as in the main text
def forsub(L,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in range(n):
        xs[i] = (bs[i] - L[i,:i]@xs[:i])/L[i,i]
    return xs

def backsub(U,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in reversed(range(n)):
        xs[i] = (bs[i] - U[i,i+1:]@xs[i+1:])/U[i,i]
    return xs

# The next two functions are the actual solution to the problem
# There is now a loop over j, since we're no longer using @
def forsub2(L,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in range(n):
        val = 0.
        for j in range(i):
            val += L[i,j]*xs[j]
        xs[i] = (bs[i] - val)/L[i,i]
    return xs

def backsub2(U,bs):
    n = bs.size
    xs = np.zeros(n)
    for i in reversed(range(n)):
        val = 0.
        for j in range(i+1,n):
            val += U[i,j]*xs[j]
        xs[i] = (bs[i] - val)/U[i,i]
    return xs

# next two functions same as in main text
def testcreate(n,val):
    A = np.arange(val,val+n*n).reshape(n,n)
    A = np.sqrt(A)
    bs = (A[0,:])**2.1
    return A, bs

def testsolve(f,A,bs):
    xs = f(A,bs); print(xs)

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    L = np.tril(A)
    testsolve(forsub,L,bs)
    testsolve(forsub2,L,bs)
    print(" ")
    U = np.triu(A)
    testsolve(backsub,U,bs)
    testsolve(backsub2,U,bs)
