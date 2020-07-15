# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from eloc import psi, ekin, epot
from montecarlo import stats
import numpy as np

def vmc(npart, ndim, al, oms, inseed=8735):
    Ncal, nm, th = 10**4, 100, 0.8
    np.random.seed(inseed)
    rolds = np.random.uniform(-1, 1, (npart, ndim))
    psiold = psi(al, oms, rolds)
    iacc, imeas = 0, 0
    eners = np.zeros(Ncal)

    for itot in range(nm*Ncal):
        rnews = rolds+th*np.random.uniform(-1,1,(npart, ndim))
        psinew = psi(al, oms, rnews)
        psiratio = (psinew/psiold)**2

        if psiratio >= np.random.uniform(0,1):
            rolds = np.copy(rnews)
            psiold = psinew
            iacc +=1
        if (itot%nm)==0:
            eners[imeas] = ekin(al,oms,rolds)+epot(oms,rolds)
            imeas += 1

    return iacc/(nm*Ncal), eners

if __name__ == '__main__':
    npart, ndim, al = 4, 3, 0.6
    oms = np.arange(1, 1 + ndim)
    accrate, eners = vmc(npart, ndim, al, oms)
    av, err = stats(eners)
    print(accrate, av, err)
