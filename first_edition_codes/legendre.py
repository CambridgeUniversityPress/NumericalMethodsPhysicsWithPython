# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

import matplotlib.pyplot as plt

def legendre(n,x):
    if n==0:
        val2 = 1.
        dval2 = 0.
    elif n==1:
        val2 = x
        dval2 = 1.
    else:
        val0 = 1.; val1 = x
        for j in range(1,n):
            val2 = ((2*j+1)*x*val1 - j*val0)/(j+1)
            val0, val1 = val1, val2
        dval2 = n*(val0-x*val1)/(1.-x**2)
    return val2, dval2

def plotlegendre(der,nsteps):
    plt.xlabel('$x$', fontsize=20)

    dertostr = {0: "$P_n(x)$", 1: "$P_n'(x)$"}
    plt.ylabel(dertostr[der], fontsize=20)
        
    ntomarker = {1: 'k-', 2: 'r--', 3: 'b-.', 4: 'g:', 5: 'c^'}
    xs = [i/nsteps for i in range (-nsteps+1,nsteps)]
    for n,marker in ntomarker.items():
        ys = [legendre(n,x)[der] for x in xs]
        labstr = 'n={0}'.format(n)
        plt.plot(xs, ys, marker, label=labstr, linewidth=3)

    plt.ylim(-3*der-1, 3*der+1)
    plt.legend(loc=4)
    plt.show()

if __name__ == '__main__':
    nsteps = 200
    plotlegendre(0,nsteps)
    plotlegendre(1,nsteps)
