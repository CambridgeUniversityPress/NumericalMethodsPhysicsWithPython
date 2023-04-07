from jacobi import termcrit
from newnormal import generatedata
from numpy import exp
import numpy as np

def soft(der, x):
    return np.log(1+exp(x)) if der==0 else exp(x)/(1+exp(x))

def getzspofx(W, x, act):
    zs = act(0, W[1,:]*x + W[2,:])
    pofx = act(0, W[0,:]@zs + W[3,0])
    return zs, pofx

def backprop(W, x, y, act):
    zs, pofx = getzspofx(W, x, act)
    ders = np.zeros((4,zs.size))
    ders[3,0] = 2*(pofx - y)*act(1, W[0,:]@zs + W[3,0])
    ders[0,:] = ders[3,0]*zs
    ders[1,:] = ders[3,0]*W[0,:]*x*act(1,W[1,:]*x+W[2,:])
    ders[2,:] = ders[3,0]*W[0,:]*act(1,W[1,:]*x+W[2,:])
    return ders

def sgd(W, x, y, act=soft, g=0.006, kmax=300, tol=1.e-4):
    for k in range(1,kmax):
        ders = backprop(W, x, y, act)
        Wn = np.array([W[i,:]-g*ders[i,:] for i in range(4)])
        _, pofxnew = getzspofx(Wn, x, act)
        conds = [termcrit(W[i,:],Wn[i,:]) for i in range(3)]
        if np.all([val<tol for val in conds[:3]]):
            break
        W = np.copy(Wn)
    return Wn, (pofxnew - y)**2

if __name__ == '__main__':
    sts = (1.5,1.2,0.15,0.08)
    data = generatedata(400,0.2,4.2,sts)
    np.random.seed(5179)
    W = np.random.uniform(-1,1,(4,10))
    W, s = sgd(W, data[0,155], data[1,155]); print(s)
