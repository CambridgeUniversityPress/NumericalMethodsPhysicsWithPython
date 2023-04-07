# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from barycentric import f, generatedata
from gauelim_pivot import gauelim_pivot
import numpy as np

def computecs(dataxs,datays):
    n = dataxs.size
    A = np.zeros((n-2,n-2))
    np.fill_diagonal(A, 2*(dataxs[2:]-dataxs[:-2]))
    np.fill_diagonal(A[1:,:], dataxs[2:-1]-dataxs[1:-2])
    np.fill_diagonal(A[:,1:], dataxs[2:-1]-dataxs[1:-2])

    b1 = (datays[2:]-datays[1:-1])/(dataxs[2:]-dataxs[1:-1])
    b2 = (datays[1:-1]-datays[:-2])/(dataxs[1:-1]-dataxs[:-2])
    bs = 6*(b1 - b2)

    cs = np.zeros(n)
    cs[1:-1] = gauelim_pivot(A, bs)
    return cs

def splineinterp(dataxs,datays,cs,x):
    k = np.argmax(dataxs>x)
    xk = dataxs[k]; xk1 = dataxs[k-1]
    yk = datays[k]; yk1 = datays[k-1]
    ck = cs[k]; ck1 = cs[k-1]

    val = yk1*(xk-x)/(xk-xk1) + yk*(x-xk1)/(xk-xk1)
    val -= ck1*((xk-x)*(xk-xk1) - (xk-x)**3/(xk-xk1))/6
    val -= ck*((x-xk1)*(xk-xk1) - (x-xk1)**3/(xk-xk1))/6
    return val

if __name__ == '__main__':
    dataxs, datays = generatedata(15, f, "equi")
    cs = computecs(dataxs, datays)
    x = 0.95; pofx = splineinterp(dataxs, datays, cs, x) 
    print(x, pofx, f(x))
