# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 3, problem 10

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

def phi(xs):
    x0, x1 = xs
    return x0**2 - 2*x0 + x1**4 - 2*x1**2 + x1

# This is chapter 3, so we are using lists (instead of numpy arrays)
def gradient(phi,xs,h=1.e-6):
    n = len(xs)
    grad = [0. for i in range(n)]
    xphs = xs[:]
    phi_xs = phi(xs)
    for j,x in enumerate(xs):
        xphs[j] = x + h
        phi_xphs = phi(xphs)
        grad[j] = (phi_xphs - phi_xs)/h
        xphs[j] = x
    return grad

npoints = 200
xs = [-0.501 + i*2.501/(npoints-1) for i in range(npoints)]
ys = [-1.7 + i*3.2/(npoints-1) for i in range(npoints)]
Z = [[0. for k in range(npoints)] for j in range(npoints)]

for i,x in enumerate(xs):
    for j,y in enumerate(ys):
        # j, i here is super important
        # https://eli.thegreenplace.net/2014/meshgrids-and-disambiguating-rows-and-columns-from-cartesian-coordinates/
        dum = gradient(phi, [x,y])
        Z[j][i] = dum[0]*dum[1]


# There are several other ways of visualizing this, which
# you may wish to experiment with.
ax = plt.axes(projection='3d')
ax.contour3D(xs, ys, Z, 80, cmap='binary')
ax.set_xlabel('$x_0$', fontsize=18)
ax.set_ylabel('$x_1$', fontsize=18)
ax.set_zlabel(r"$(\partial \phi/\partial x)(\partial \phi/\partial y)$", fontsize=18)

ax.view_init(34, -168)

plt.show()
