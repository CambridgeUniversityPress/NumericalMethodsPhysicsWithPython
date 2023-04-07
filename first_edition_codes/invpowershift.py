# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triang import forsub, backsub, testcreate
from ludec import ludec
from jacobi import termcrit
from power import mag, testeigone
import numpy as np

#def invpowershift(A,shift=20,kmax=200,tol=1.e-2):
def invpowershift(A,shift=20,kmax=200,tol=1.e-8):
    n = A.shape[0]
    znews = np.ones(n)
    qnews = znews/mag(znews)
    Astar = A - np.identity(n)*shift
    L, U = ludec(Astar)

    for k in range(1,kmax):
        qs = np.copy(qnews)
        ys = forsub(L,qs)
        znews = backsub(U,ys)
        qnews = znews/mag(znews)

        if qs@qnews<0:
            qnews = -qnews

        err = termcrit(qs,qnews)
        print(k, qnews, err)

        if err < tol:
            lam = qnews@A@qnews
            break
    else:
        lam = qnews = None

    return lam, qnews

if __name__ == '__main__':
    A, bs = testcreate(4,21)
    testeigone(invpowershift,A)
