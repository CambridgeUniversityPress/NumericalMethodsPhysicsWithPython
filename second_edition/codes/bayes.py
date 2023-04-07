from gauelim_pivot import gauelim_pivot
from newnormal import generatedata, phi
import numpy as np

def inv(A):
    n = A.shape[0]
    invA = np.zeros((n,n))
    for i,bs in enumerate(np.identity(n).T):
        invA[:,i] = gauelim_pivot(A,bs)
    return invA

def bayes(data, batches, primus, priS):
    n = primus.size
    i = 0
    for N in batches:
        A = np.zeros((N,n))
        for k in range(n):
            A[:,k] = phi(n,k,data[0,i:i+N])/data[2,i:i+N]
        bs = data[1,i:i+N]/data[2,i:i+N]

        priSinv = inv(priS)
        postS = inv(A.T@A + priSinv)
        postmus = postS@(A.T@bs + priSinv@primus)
        
        primus, priS = postmus, postS
        i += N
        print(n, i, postmus, postS)
    return postmus, postS
    
if __name__ == '__main__':
    batches = 1, 4, 3
    data = generatedata(np.sum(batches))
    for n in (5, 2):
        primus = np.zeros(n)
        priS = np.zeros((n,n))
        np.fill_diagonal(priS, np.linspace(10,2,n))
        postmus, postS = bayes(data, batches, primus, priS)
        print(postmus)
