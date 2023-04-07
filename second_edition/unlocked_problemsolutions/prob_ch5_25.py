# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (2nd ed., CUP, 2023)

# Solution to chapter 5, problem 25

# NOTE TO INSTRUCTORS: this solution is made available to all readers (i.e., not locked)

from math import sqrt, exp

def f(x):
    # following line is there for development/debugging purposes
    # that is, so you can explicitly check when the function is called.
    print("FUNCTION EVALUATION, x= ", x)
    return exp(x - sqrt(x)) - x


# this is similar to golden() but we save intermediate variables
# Also, depending on the conditional/if branch selected, we 
# store different things. As always, the more efficient version
# of the code is messier/uglier
def golden(f,x0,x1,kmax=200,tol=1.e-8):
    phi = 0.5*(1 + sqrt(5))
    x2 = x0 + (x1-x0)/(phi+1)
    f2 = f(x2)
    x3 = x0 + (x1-x0)*phi/(phi+1)
    f3 = f(x3)
    for k in range(1,kmax):
        if f2<f3:
            x1 = x3
            x3,f3 = x2,f2
            x2 = x0 + (x1-x0)/(phi+1)
            f2 = f(x2)
        else:
            x0 = x2
            x2,f2 = x3,f3
            x3 = x0 + (x1-x0)*phi/(phi+1)
            f3 = f(x3)
            
        xnew = (x0+x1)/2
        xdiff = abs(x1-x0)
        rowf = "{0:2d} {1:1.16f} {2:1.16f} {3:1.16f}"
        #print(rowf.format(k,xnew,xdiff,f(xnew)))

        if abs(xdiff) < tol:
            break
    else:
        xnew = None

    return xnew

if __name__ == '__main__':
    val = golden(f,0.,3.5)
    print(val)
