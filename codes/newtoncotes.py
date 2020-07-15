# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import numpy as np

def f(x):
    return 1/np.sqrt(x**2 + 1)

def rectangle(f,a,b,n):
    h = (b-a)/(n-1)
    xs = a + np.arange(n-1)*h
    fs = f(xs)
    return h*np.sum(fs)

def trapezoid(f,a,b,n):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    cs = np.ones(n); cs[0] = 0.5; cs[-1] = 0.5
    contribs = cs*f(xs)
    return h*np.sum(contribs)

def simpson(f,a,b,n):
    h = (b-a)/(n-1)
    xs = a + np.arange(n)*h
    cs = 2*np.ones(n)
    cs[1::2] = 4; cs[0] = 1; cs[-1] = 1
    contribs = cs*f(xs)
    return (h/3)*np.sum(contribs)

if __name__ == '__main__':
    ans = np.log(1 + np.sqrt(2))
    print(ans)

    for integrator in (rectangle, trapezoid, simpson):
        print(integrator(f, 0., 1., 51), end=" ")
