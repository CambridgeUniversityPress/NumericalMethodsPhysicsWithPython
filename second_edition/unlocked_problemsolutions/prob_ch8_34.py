# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 8, problem 34

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

# this is very similar to matsetup in evp_matrix.py
def Hmat_order2(n,xmax):
    a = 2*xmax/(n-1)
    val = 0.5/a**2
    xs = np.linspace(-xmax,+xmax,n)
    #print(xs)

    Hmat = np.zeros((n,n))
    np.fill_diagonal(Hmat, 0.5*xs**2)
    Hmat += 2*val*np.eye(n)
    Hmat -= val*np.eye(n,k=1)
    Hmat -= val*np.eye(n,k=-1)
    return Hmat

if __name__ == '__main__':

    # We use NumPy's eigenvalues
    n = 50; xmax = 10.
    Hmat = Hmat_order2(n,xmax)
    eigvalsA, eigvecs = eig(Hmat)

    n = 100; xmax = 10.
    Hmat = Hmat_order2(n,xmax)
    eigvalsB, eigvecs = eig(Hmat)

    # The only subtlety is in figuring out the order
    # in which the eigenvalues are stored.
    # To avoid getting lost, we simply sort the eigenvalues:
    eigA = np.sort(eigvalsA)
    eigB = np.sort(eigvalsB)
    for i in range(6):
        print(eigB[i])
        extr = (eigB[i]*4 - eigA[i])/3
        print(extr)
        print("")
    # The extrapolated values are all correct to three significant figures.
    # This is a considerable improvement, even though we used the same discretization (twice).
