# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from kron import paulimatrices, kron
from qrmet import qrmet
import numpy as np

def twospins(omI,omII,gam):
    hbar = 1.
    paulis = paulimatrices()
    iden = np.identity(2)

    SIs = [hbar*kron(pa,iden)/2 for pa in paulis]
    SIIs = [hbar*kron(iden,pa)/2 for pa in paulis]
    SIdotII = sum([SIs[i]@SIIs[i] for i in range(3)])

    H = -omI*SIs[2] - omII*SIIs[2] + gam*SIdotII
    H = H.real
    return H

if __name__ == '__main__':
    H = twospins(1.,2.,0.5)
    qreigvals = qrmet(H); print(qreigvals)
