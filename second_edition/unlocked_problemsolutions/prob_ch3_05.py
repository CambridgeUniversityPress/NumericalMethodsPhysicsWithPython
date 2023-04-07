# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 3, problem 5

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import exp, sin, cos, log10
import matplotlib.pyplot as plt

def f(x):
    return exp(sin(2*x))

def fprime(x):
    return 2*exp(sin(2*x))*cos(2*x)

def calc_fd(f,x,h):
    fd = (f(x+h) - f(x))/h
    return fd

def calc_cd(f,x,h):
    cd = (f(x+h/2) - f(x-h/2))/h
    return cd

def calc_fd2(f,x,h):
    fd2 = (4*f(x+h/2) - f(x+h) - 3*f(x))/h
    return fd2

def calc_cd2(f,x,h):
    cd2 = (27*f(x+h/2) + f(x-3*h/2) - 27*f(x-h/2) - f(x+3*h/2))/(24*h)
    return cd2

def logplot(hs,fds,cds,fds2,cds2):
    plt.xlabel('log(h)')
    plt.ylabel('log(|abs. error|)')

    xs = [log10(h) for h in hs]
    y1s = [log10(fd) for fd in fds]
    y2s = [log10(cd) for cd in cds]
    y3s = [log10(fd2) for fd2 in fds2]
    y4s = [log10(cd2) for cd2 in cds2]

    plt.plot(xs,y1s,'r-o',label='forward diff.')
    plt.plot(xs,y2s,'b--s',label='central diff.')
    plt.plot(xs,y3s,'g-.*',label='forward diff. 2')
    plt.plot(xs,y4s,'y:D',label='central diff. 2')

    plt.legend()
    plt.show()

if __name__ == '__main__':
    x = 0.5
    an = fprime(x)

    hs = [10**(-i) for i in range(1,18)]
    fds = [abs(calc_fd(f,x,h) - an) for h in hs]
    cds = [abs(calc_cd(f,x,h) - an) for h in hs]
    fds2 = [abs(calc_fd2(f,x,h) - an) for h in hs]
    cds2 = [abs(calc_cd2(f,x,h) - an) for h in hs]

    rowformat = "{0:1.0e} {1:1.16f} {2:1.16f}"
    print("h     abs. error in fd   abs. error in cd")
    for h,fd,cd in zip(hs,fds,cds):
        print(rowformat.format(h,fd,cd))

    print("h     abs. error in fd2  abs. error in cd2")
    for h,fd2,cd2 in zip(hs,fds2,cds2):
        print(rowformat.format(h,fd2,cd2))

    logplot(hs,fds,cds,fds2,cds2)

#The second forward difference is basically as good as the first central difference. 
#This was to be expected, based on the approximation error dependence: O(h^2)
#The second central difference is way better than everything else: it gives
#a larger h_{opt} \approx 10^{-3} that leads to a minimum absolute error that is 
#smaller than the other ones by a couple of orders of magnitude.
#The second forward difference improves by 2 orders for each H step and the 
#second central difference improves by 4 orders, before roundoff starts dominating
#and then it does as poorly/well as everybody else.
