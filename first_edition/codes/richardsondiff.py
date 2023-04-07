# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

from finitediff import f, fprime, calc_fd, calc_cd

x = 0.5
an = fprime(x)

hs = [10**(-i) for i in range(1,7)]

rowf = "{0:1.0e} {1:1.16f} {2:1.16f}"
print("h     abs. err. rich fd  abs. err. rich cd")
for h in hs:
    fdrich = 2*calc_fd(f,x,h/2) - calc_fd(f,x,h)
    fd = abs(fdrich-an)
    cdrich = (4*calc_cd(f,x,h/2) - calc_cd(f,x,h))/3
    cd = abs(cdrich-an)
    print(rowf.format(h,fd,cd))
