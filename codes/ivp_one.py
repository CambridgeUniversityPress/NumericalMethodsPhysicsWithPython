# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def f(x,y):
    return - (30/(1-x**2)) + ((2*x)/(1-x**2))*y - y**2

def euler(f,a,b,n,yinit):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros(n)

    y = yinit
    for j,x in enumerate(xs):
        ys[j] = y
        y += h*f(x, y)
    return xs, ys
        
def rk4(f,a,b,n,yinit):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    ys = np.zeros(n)

    y = yinit
    for j,x in enumerate(xs):
        ys[j] = y
        k0 = h*f(x, y)
        k1 = h*f(x+h/2, y+k0/2)
        k2 = h*f(x+h/2, y+k1/2)
        k3 = h*f(x+h, y+k2)
        y += (k0 + 2*k1 + 2*k2 + k3)/6
    return xs, ys
        
if __name__ == '__main__':
    a, b, n, yinit = 0.05, 0.49, 12, 19.53
    xs, ys = euler(f,a,b,n,yinit); print(ys)
    xs, ys = rk4(f,a,b,n,yinit); print(ys)
