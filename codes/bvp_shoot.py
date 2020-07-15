# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from secant import secant
from ivp_two import fs, rk4_gen
import numpy as np

def shoot(sig):
    a, b, n = 0.05, 0.49, 100
    yinits = np.array([0.0926587109375, sig])
    xs, ws = rk4_gen(fs,a,b,n,yinits)  
    wfinal = 0.11177050858750004
    return ws[-1, 0] - wfinal

if __name__ == '__main__':
    wder = secant(shoot,0.,1.) 
    print(wder)
