# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from kahansum import kahansum
from math import sqrt

def chargearray(nvals):
    vals = [-0.5 + i/(nvals-1) for i in range(nvals)]
    qtopos = {}
    for i,posx in enumerate(vals):
        for j,posy in enumerate(vals):
            count = j + nvals*i + 1
            key = 1.02*count if (i+j)%2==0 else -count
            qtopos[key] = posx, posy
    return qtopos

def vecmag(rs):
    sq = [r**2 for r in rs]
    return sqrt(kahansum(sq))

def fullpot(qtopos,rs):
    potvals = []
    for q,pos in qtopos.items():
        diffs = [r - po for r,po in zip(rs,pos)]
        R = vecmag(diffs)
        potvals.append(q/R)
    return kahansum(potvals)

if __name__ == '__main__':
    qtopos = chargearray(6)
    for y in 1,-1:
        rs = [0.,y]
        potval = fullpot(qtopos,rs) 
        print(rs, potval)
