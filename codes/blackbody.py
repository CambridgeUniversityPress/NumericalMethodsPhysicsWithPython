# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from multi_newton import jacobian, multi_newton
import numpy as np

def generatedata():
    dataxs = np.array([373.1, 492.5, 733, 755, 799, 820,
            877, 1106, 1125, 1403, 1492, 1522, 1561])
    datays = np.array([156., 638, 3320, 3810, 4440, 5150,
            6910, 16400, 17700, 44700, 57400, 60600, 67800])
    datasigs = np.ones(dataxs.size)
    return dataxs, datays, datasigs

def model(cs,x):
    return cs[0] + cs[1]*x**cs[2]

def fs(cs):
    dataxs, datays, datasigs = generatedata()
    c0, c1, c2 = cs
    resids = datays - model(cs, dataxs)
    f0 = np.sum(resids/datasigs**2)
    f1 = np.sum(dataxs**c2*resids/datasigs**2)
    numers = c1*dataxs**c2*np.log(dataxs)*resids
    f2 = np.sum(numers/datasigs**2)
    return np.array([f0,f1,f2])

def computechisq(cs):
    dataxs, datays, datasigs = generatedata()
    chisq = np.sum((datays - model(cs,dataxs))**2/datasigs**2)
    return chisq

if __name__ == '__main__':
    colds = np.array([-700, 1.26e-8, 6])
    cs = multi_newton(fs,jacobian,colds,kmax=500,tol=2.e-6)
    chisq = computechisq(cs); print(cs); print(chisq)
