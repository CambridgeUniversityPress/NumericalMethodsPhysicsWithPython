# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from newtoncotes import f, trapezoid
import numpy as np

def prettyprint(row):
    for elem in row:
        print("{0:1.10f} ".format(elem),end="")
    print("")

def richardson(Rprev, Rincurr0, i):
    Rcurr = [Rincurr0]
    for j in range(1, i+1):
        val = (4**j*Rcurr[j-1] - Rprev[j-1])/(4**j - 1)
        Rcurr.append(val)
    return Rcurr

def romberg(f,a,b,imax = 20,tol = 1.e-8):
    n = 2
    val = trapezoid(f,a,b,n) 
    Rprev = [val]
    #prettyprint(Rprev)
    for i in range(1,imax):
        nprime = 2*n-1
        Rincurr0 = trapezoid(f,a,b,nprime)
        Rcurr = richardson(Rprev, Rincurr0, i)
        #prettyprint(Rcurr)
        err = abs(Rprev[-1] - Rcurr[-1])/abs(Rcurr[-1])
        valprime = Rcurr[-1]
        if err < tol:
            break
        n = nprime
        Rprev = Rcurr[:]
    else:
        valprime = None
    return valprime

if __name__ == '__main__':
    ans = np.log(1 + np.sqrt(2))
    print(ans)
    print(romberg(f,0.,1.))
