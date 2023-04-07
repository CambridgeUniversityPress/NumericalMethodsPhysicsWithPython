# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

from eloc import psi, ekin, epot
from montecarlo import stats
import numpy as np

def vmc(npart, ndim, al, oms, inseed=8735):
    ncurly, nm, th = 10**4, 100, 0.8
    np.random.seed(inseed)
    Xold = np.random.uniform(-1, 1, (npart, ndim))
    psiold = psi(al, oms, Xold)
    iacc, imeas = 0, 0
    eners = np.zeros(ncurly)

    for itot in range(nm*ncurly):
        Xnew = Xold + th*np.random.uniform(-1,1,(npart, ndim))
        psinew = psi(al, oms, Xnew)
        psiratio = (psinew/psiold)**2

        if np.random.uniform(0,1) <= psiratio:
            Xold = np.copy(Xnew)
            psiold = psinew
            iacc +=1
        if (itot%nm)==0:
            eners[imeas] = ekin(al,oms,Xold) + epot(oms,Xold)
            imeas += 1

    return iacc/(nm*ncurly), eners

if __name__ == '__main__':
    npart, ndim, al = 4, 3, 0.6
    oms = np.arange(1, 1 + ndim)
    accrate, eners = vmc(npart, ndim, al, oms)
    av, err = stats(eners)
    print(accrate, av, err)
