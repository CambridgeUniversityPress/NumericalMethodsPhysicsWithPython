# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from qrmet import qrmet
from kron import paulimatrices, kron
import numpy as np

def threespins(omI,omII,omIII,gam):
    hbar = 1.
    paulis = paulimatrices()
    iden = np.identity(2)

    SIs = [hbar*kron(kron(pa,iden),iden)/2 for pa in paulis]
    SIIs = [hbar*kron(kron(iden,pa),iden)/2 for pa in paulis]
    SIIIs = [hbar*kron(kron(iden,iden),pa)/2 for pa in paulis]

    SIdotII = sum([SIs[i]@SIIs[i] for i in range(3)])
    SIdotIII = sum([SIs[i]@SIIIs[i] for i in range(3)])
    SIIdotIII = sum([SIIs[i]@SIIIs[i] for i in range(3)])

    H = -omI*SIs[2] - omII*SIIs[2] - omIII*SIIIs[2]
    H += gam*(SIdotII+SIdotIII+SIIdotIII)
    H = H.real
    return H

if __name__ == '__main__':
    np.set_printoptions(precision=3)
    H = threespins(1.,2.,3.,0.5)
    qreigvals = qrmet(H); print(qreigvals)
