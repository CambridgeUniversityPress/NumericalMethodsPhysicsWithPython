# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 4, problem 49

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

# We import a lot of functionality from the text codes
import numpy as np
from power import mag
from triang import forsub, backsub
from ludec import ludec
from eig import testeigall
from qrmet import qrmet
from kron import paulimatrices, kron
from twospins import twospins
from threespins import threespins

# As desired, this does not call termcrit(), but 
# simply carries out the same operation over and over 
# without checking for convergence.
# Step is necessary because err calculation in invpowershift() is misleading: 
# eigvectors are good, but when you are dividing with tiny numbers it looks like it's large.
# Conceptually, this is because the eigenvectors become (1,0,0,0) etc, i.e., they have many zeros
def invpowershift_notfancy(A,shift=20,kmax=200):
    n = A.shape[0]
    znews = np.ones(n)
    qnews = znews/mag(znews)
    Astar = A - np.identity(n)*shift
    L, U = ludec(Astar)

    for k in range(kmax):
        qs = np.copy(qnews)
        ys = forsub(L,qs)
        znews = backsub(U,ys)
        qnews = znews/mag(znews)
        #print(k,qnews)

    Aqnews = np.dot(A,qnews) 
    lam = np.dot(qnews,Aqnews)
    return lam, qnews

# This is identical to the eig() from eig.py, but it
# calls our new inverse-power method above.
def eigmod(A,eps=1.e-12):
    n = A.shape[0]
    eigvals = np.zeros(n)
    eigvecs = np.zeros((n,n))
    qreigvals = qrmet(A)
    for i, qreigval in enumerate(qreigvals):
        eigvals[i], eigvecs[:,i] = invpowershift_notfancy(A,qreigval+eps)
    return eigvals, eigvecs

if __name__ == '__main__':
    np.set_printoptions(precision=3)
    # We are benchmarking the eigenvectors against NumPy

    A = twospins(1.,2.,0.5)
    eigvals, eigvecs = eigmod(A)
    print(eigvals)
    print(eigvecs)
    print("")
    npeigvals, npeigvecs = np.linalg.eig(A)
    print(npeigvals)
    print(npeigvecs)
    print("")
    print("")

    A = threespins(1.,2.,3.,0.5)
    eigvals, eigvecs = eigmod(A)
    print(eigvals)
    print(eigvecs)
    print("")
    npeigvals, npeigvecs = np.linalg.eig(A)
    print(npeigvals)
    print(npeigvecs)
