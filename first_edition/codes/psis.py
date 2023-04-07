# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from math import sqrt, pi, factorial, exp
import cmath

def hermite(n,x):
    val0 = 1.; val1 = 2*x
    for j in range(1,n):
        val2 = 2*x*val1 - 2*j*val0
        val0, val1 = val1, val2
    dval2 = 2*n*val0
    return val2, dval2

def psiqho(x,nametoval):
    n = nametoval["n"]
    momohbar = nametoval["momohbar"]
    al = nametoval["al"]
    psival = momohbar**0.25*exp(-0.5*al*momohbar * x**2) 
    psival *= hermite(n,sqrt(momohbar)*x)[0]
    psival /= sqrt(2**n * factorial(n) * sqrt(pi))
    return psival

def psibox(x,nametoval):
    n = nametoval["n"]
    boxl = nametoval["boxl"]
    return cmath.exp(2*pi*n*x*1j/boxl)

if __name__ == '__main__':
    x = 1.
    nametoval = {"n": 100, "momohbar": 1., "al": 1.}
    psiA = psiqho(x, nametoval)
    nametoval = {"n": -2, "boxl": 2*pi}
    psiB = psibox(x, nametoval)
    print(psiA, psiB)
