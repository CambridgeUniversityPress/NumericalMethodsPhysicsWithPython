# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from fft import fft
import numpy as np

def fft2(A):
    B = A.astype(complex)
    for i, row in enumerate(A):
        B[i,:] = fft(row)
    for j, col in enumerate(B.T):
        B[:,j] = fft(col)
    return B

def inversefft2(A):
    n2 = A.size
    newA = fft2(np.conjugate(A))/n2
    return np.conjugate(newA)

def func(x,y):
    return np.cos(3*x+4*y) - np.cos(5*x-2*y)

def poisson(a,b,n):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    Xs, Ys = np.meshgrid(xs,xs)

    F = func(Xs, Ys)
    Ftil = fft2(F)

    ks = np.arange(n)
    Kxs, Kys = np.meshgrid(ks,ks)
    Denom = np.cos(2*np.pi*Kxs/n) + np.cos(2*np.pi*Kys/n) - 2
    Phitil = 0.5*Ftil*h**2/Denom
    Phitil[0,0] = 0

    Phi = np.real(inversefft2(Phitil))
    return Phi

if __name__ == '__main__':
    Phi = poisson(0, 2*np.pi, 128); print(Phi)
