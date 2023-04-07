from triang import testcreate
from invpowershift import invpowershift
from qrmet import qrmet
import numpy as np

def makeH(A):
    n = A.shape[0]
    H = np.zeros((2*n,2*n))
    H[n:,:n] = A
    H[:n,n:] = A.T
    return H

def svd(A, solver=qrmet):
    H = makeH(A)
    shift = 1.
    Hstar = H - np.identity(H.shape[0])*shift
    Hvals = solver(Hstar) + shift
    indices = (Hvals>=0).nonzero()[0]
    S = Hvals[indices]

    n = A.shape[0]
    vals = np.zeros(2*n)
    vecs = np.zeros((2*n,2*n))
    for i, qre in enumerate(Hvals):
        vals[i], vecs[:,i] = invpowershift(H,qre+1.e-20,tol=1.e-8)
    Vs = [vecs[:n,i] for i in indices]
    V = np.sqrt(2)*np.column_stack(Vs)
    Us = [vecs[n:,i] for i in indices]
    U = np.sqrt(2)*np.column_stack(Us)
    return U, np.diag(S), V.T 

if __name__ == '__main__':
    A, _ = testcreate(4,21)
    A -= np.identity(4)
    U, S, VT = svd(A)
    print(np.all(np.linalg.norm(A - U@S@VT)<1.e-6), S)
    npU, npS, npVT = np.linalg.svd(A)
    diffA =  A - npU@np.diag(npS)@npVT
    print(np.all(np.linalg.norm(diffA)<1.e-6), npS)
