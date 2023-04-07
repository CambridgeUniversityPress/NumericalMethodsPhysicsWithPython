# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from gauelim_pivot import gauelim_pivot
from jacobi import termcrit
import numpy as np

def generatedata():
    data = np.zeros((3,13))
    data[0,:] = np.array([373.1, 492.5, 733, 755, 799, 820,
            877, 1106, 1125, 1403, 1492, 1522, 1561])
    data[1,:] = np.array([156., 638, 3320, 3810, 4440, 5150,
            6910, 16400, 17700, 44700, 57400, 60600, 67800])
    data[2,:] = np.ones(data.shape[1])
    return data

def model(cs,x):
    return cs[0] + cs[1]*x**cs[2]

def getKrs(data, cs):
    K = np.zeros((data.shape[1], cs.size))
    K[:,0] = 1/data[2,:]
    K[:,1] = data[0,:]**cs[2]/data[2,:]
    K[:,2] = cs[1]*data[0,:]**cs[2]*np.log(data[0,:])/data[2,:]
    rs = (data[1,:] - model(cs, data[0,:]))/data[2,:]
    return K, rs

def gaussnewton(data, colds, kmax=50, tol=1.e-8):
    for k in range(1,kmax):
        K, rs = getKrs(data, colds)
        cnews = colds + gauelim_pivot(K.T@K, K.T@rs)
        err = termcrit(colds,cnews)
        print(k,err)
        if err < tol:
            break
        colds = np.copy(cnews)
    return cnews

if __name__ == '__main__':
    data = generatedata()
    colds = np.array([-700, 1.26e-8, 4.5])
    cs = gaussnewton(data,colds); print(cs)
