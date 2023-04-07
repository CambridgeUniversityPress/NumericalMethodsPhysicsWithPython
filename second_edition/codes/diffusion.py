from triang import forsub, backsub
from ludec import ludec
import numpy as np

def ftcs(fac, m, polds):
    pnews = np.copy(polds)
    for k in range(m-1):
        pnews[1:-1] = fac*(polds[:-2] + polds[2:] - 2*polds[1:-1])
        pnews[1:-1] += polds[1:-1]
        polds = np.copy(pnews)
    return pnews

def btcs(fac, m, polds):
    n = polds.size
    A = np.zeros((n,n))
    np.fill_diagonal(A, 1.+2*fac)
    A[0,0] = 1; A[-1,-1] = 1
    np.fill_diagonal(A[1:,:], -fac)   
    A[-1,-2] = 0
    np.fill_diagonal(A[:,1:], -fac)
    A[0,1] = 0
    LA, UA = ludec(A)

    for k in range(m-1):
        ys = forsub(LA,polds)
        pnews = backsub(UA,ys)
        polds = np.copy(pnews)

    return pnews

if __name__ == '__main__':
    al, L, T, n, m = 0.1, 1., 0.4, 51, 200
    hx = L/(n-1); ht = T/(m-1)
    fac = al*ht/hx**2
    pinits = 30*np.ones(n)
    pinits[:n//2] = 10
    phis = ftcs(fac, m, pinits); print(phis)
    phis = btcs(fac, m, pinits); print(phis)
