# Author: Alex Gezerlis
# Numerical Methods in Physics with Python (CUP, 2020)

def kahansum(xs):
    s = 0.; e = 0.
    for x in xs:
        temp = s 
        y = x + e
        s = temp + y
        e = (temp - s) + y
    return s
    
if __name__ == '__main__':
    xs = [0.7, 0.1, 0.3]
    print(sum(xs), kahansum(xs))
