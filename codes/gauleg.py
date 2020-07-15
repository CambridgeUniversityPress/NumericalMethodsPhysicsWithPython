# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from legendre import legendre
from legroots import legroots
from newtoncotes import f
import numpy as np

def gauleg_params(n):
    xs = legroots(n)
    cs = 2/((1-xs**2)*legendre(n,xs)[1]**2)
    return xs, cs

def gauleg(f,a,b,n):
    xs, cs = gauleg_params(n)
    coeffp = 0.5*(b+a) 
    coeffm = 0.5*(b-a)
    ts = coeffp + coeffm*xs
    contribs = cs*f(ts)
    return coeffm*np.sum(contribs)

if __name__ == '__main__':
    ans = np.log(1 + np.sqrt(2))
    print(ans)
    for n in range(2,10):
        print(n, gauleg(f,0.,1,n))
