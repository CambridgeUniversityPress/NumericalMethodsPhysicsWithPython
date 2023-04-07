# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 6, problem 46

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from gauleg import gauleg_params
import numpy as np
from scipy.stats import chi2

# The only subtlety in the following functions is that we have
# chosen to pass in nu as an argument and therefore modified
# our program from gauleg.py. If you want, instead, to use
# np.polynomial.legendre.leggauss(), which doesn't accept
# extra arguments, you could simply hardcode a given nu value in
# new function definition (via def etc.)
# or employ a lambda that corresponds to a given nu value

def parta(y, nu):
    val = np.exp(-y/2) * y**(nu/2)
    val = val/(2**(nu/2)*np.math.factorial((nu/2)-1))
    return val

def partb(y, nu):
    val = np.exp(-y/2) * y**(nu/2 - 1)
    val *= (y-nu)**2
    val = val/(2**(nu/2)*np.math.factorial((nu/2)-1))
    return val

def partc_integrand(y, nu):
    val = np.exp(-y/2) * y**(nu/2 - 1)
    val = val/(2**(nu/2)*np.math.factorial((nu/2)-1))
    return val

def gauleg(f,a,b,n, nu):
    xs, cs = gauleg_params(n)
    coeffp = 0.5*(b+a) 
    coeffm = 0.5*(b-a)
    ts = coeffp + coeffm*xs
    contribs = cs*f(ts, nu)
    return coeffm*np.sum(contribs)

# This is called adhoc because we are passing in an extra argument (nu),
# just like in the earlier functions
def adhoc_bisection(nu,x0,x1,kmax=200,tol=1.e-8):
    f0 = gauleg(partc_integrand,x0,1000.,200, nu) - 0.5
    for k in range(1,kmax):
        x2 = (x0+x1)/2
        f2 = gauleg(partc_integrand,x2,1000.,200, nu) - 0.5
        
        if f0*f2 < 0:
            x1 = x2
        else:
            x0, f0 = x2, f2
            
        x2new = (x0+x1)/2
        xdiff = abs(x2new-x2)
        rowf = "{0:2d} {1:1.16f} {2:1.16f} {3:1.16f}"

        if abs(xdiff/x2new) < tol:
            break
    else:
        x2new = None

    return x2new

def partd(y, nu):
    val = np.exp(-y/2) * y**(nu/2 - 1)
    val = val/(2**(nu/2)*np.math.factorial((nu/2)-1))
    return val

if __name__ == '__main__':
    # Note that we are using Gauss-Legendre quadrature, which integrates over a finite interval.
    # We therefore pass in a finite (but large) right endpoint, b.
    # part (a)
    for nu in 4, 8, 12:
        print(nu, gauleg(parta,0.,1000.,200, nu))
    # part (b)
    for nu in 4, 8, 12:
        print(2*nu, gauleg(partb,0.,1000.,200, nu))

    # you may have to experiment a little on the input guesses given to the bisection function
    # part (c)
    for nu in 4, 8, 12:
        print(nu*(1 - 2/(9*nu))**3, adhoc_bisection(nu,0.1,50), chi2.median(nu))
    # part (d)
    for M in 1,4,9:
        for nu in 4, 8, 12:
            print(M, nu, gauleg(partd,0.,M,200, nu), chi2.cdf(M, nu))
