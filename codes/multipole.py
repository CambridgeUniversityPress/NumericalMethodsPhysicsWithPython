# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from kahansum import kahansum
from chargearray import chargearray, vecmag
from legendre import legendre

def decomp(rs,ris):
    rmag = vecmag(rs); rimag = vecmag(ris)
    prs = [r*ri for r,ri in zip(rs,ris)]
    vecdot = kahansum(prs)
    costheta = vecdot/(rmag*rimag)
    return rmag, rimag, costheta

def multicoes(rs,qtopos,nmax=60):
    coes = [0. for n in range(nmax+1)]
    for n in range(nmax+1):
        for q,pos in qtopos.items():
            rmag, rimag, costheta = decomp(rs,pos)
            val = q*(rimag**n)*legendre(n,costheta)[0]
            coes[n] += val
    return coes

def multifullpot(rs,qtopos):
    coes = multicoes(rs,qtopos)
    rmag = vecmag(rs)
    contribs = [coe/rmag**(n+1) for n,coe in enumerate(coes)]
    return kahansum(contribs)

if __name__ == '__main__':
    qtopos = chargearray(6)
    for y in 1,-1:
        rs = [0.,y]
        potval = multifullpot(rs,qtopos) 
        print(rs, potval)
