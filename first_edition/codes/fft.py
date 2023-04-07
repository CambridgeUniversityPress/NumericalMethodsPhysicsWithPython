# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from triginterp import f, generatedata
from math import pi
import numpy as np

def fft(ys):
    n = ys.size
    m = n//2
    if n==1:
        ytils = ys
    else:
        evens = fft(ys[::2])
        odds = fft(ys[1::2])
        coeffs = np.exp(-2*pi*np.arange(m)*1j/n)
        first = evens + coeffs*odds
        second = evens - coeffs*odds
        ytils = np.concatenate((first, second))
    return ytils

def fftinterp(ytils,x):
    n = ytils.size
    m = n//2
    val = ytils[:m]@np.exp(np.arange(m)*x*1j)
    val += ytils[m]*np.cos(m*x)
    val += ytils[m+1:]@np.exp(np.arange(-m+1,0)*x*1j)
    return val/n

if __name__ == '__main__':
    n = 8
    dataxs, datays = generatedata(n, f)
    ytils = fft(datays)
    x = 0.3; pofx = fftinterp(ytils, x)
    print(x,pofx.real)
