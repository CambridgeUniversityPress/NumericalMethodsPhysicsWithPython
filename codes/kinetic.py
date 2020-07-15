# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from psis import psiqho, psibox
from math import pi

def kinetic(psi,x,nametoval,h=0.005):
    hom = 1.
    psiold = psi(x,nametoval)
    psip = psi(x+h,nametoval)
    psim = psi(x-h,nametoval)

    lapl = (psip + psim - 2.*psiold)/h**2
    kin = -0.5*hom*lapl/psiold
    return kin

def test_kinetic():
    x = 1.
    hs = [10**(-i) for i in range(1,6)]

    nametoval = {"n": 100, "momohbar": 1., "al": 1.}
    qhos = [kinetic(psiqho,x,nametoval,h) for h in hs]
    nametoval = {"n": -2, "boxl": 2*pi}
    boxs = [kinetic(psibox,x,nametoval,h) for h in hs]

    rowf = "{0:1.0e} {1:1.16f} {2:1.16f}"
    print("h     qho                 box")
    for h,qho,box in zip(hs,qhos,boxs):
        print(rowf.format(h,qho,box))

if __name__ == '__main__':
    test_kinetic()
