# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def fs(x,yvals):
    y0, y1 = yvals
    f0 = y1
    f1 = - (30/(1-x**2))*y0 + ((2*x)/(1-x**2))*y1 
    return np.array([f0, f1])

def rk4_gen(fs,a,b,n,yinits):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros((n, yinits.size))

    yvals = np.copy(yinits)
    for j,x in enumerate(xs):
        ys[j,:] = yvals
        k0 = h*fs(x, yvals)
        k1 = h*fs(x+h/2, yvals+k0/2)
        k2 = h*fs(x+h/2, yvals+k1/2)
        k3 = h*fs(x+h, yvals+k2)
        yvals += (k0 + 2*k1 + 2*k2 + k3)/6
    return xs, ys
        
if __name__ == '__main__':
    a, b, n = 0.05, 0.49, 12
    yinits = np.array([0.0926587109375, 1.80962109375])
    xs, ys = rk4_gen(fs,a,b,n,yinits)
    print(ys)
