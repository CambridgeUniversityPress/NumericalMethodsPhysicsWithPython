# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from secant import secant
import numpy as np

def fs(x,yvals,s):
    q = 1.5
    y0, y1 = yvals
    f0 = y1
    f1 = (2*q*np.cos(2*x) - s)*y0
    return np.array([f0, f1])

def rk4_gen_eig(fs,a,b,n,yinits,s):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros((n, yinits.size))

    yvals = np.copy(yinits)
    for j,x in enumerate(xs):
        ys[j,:] = yvals
        k0 = h*fs(x, yvals,s)
        k1 = h*fs(x+h/2, yvals+k0/2,s)
        k2 = h*fs(x+h/2, yvals+k1/2,s)
        k3 = h*fs(x+h, yvals+k2,s)
        yvals += (k0 + 2*k1 + 2*k2 + k3)/6
    return xs, ys

def shoot(s):
    a, b, n = 0, 2*np.pi, 500
    yinits = np.array([0., 5.])

    xs, ys = rk4_gen_eig(fs,a,b,n,yinits,s)  
    wfinal = 0.
    return ys[-1, 0] - wfinal

if __name__ == '__main__':
    for sinit in (-0.4, 3.3, 8.5):
        sval = secant(shoot,sinit,sinit+0.5) 
        print(sval, end=" ")
