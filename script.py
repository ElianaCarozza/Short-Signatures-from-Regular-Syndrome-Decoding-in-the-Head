import math 
from scipy.misc import comb
def Func_u(N,H,x,i,j): # x corresponds to theta in the writeup , and H to ell in the writeup
    val = (6 **j)*comb(H,j,exact=True)*comb(N/6-H,i-j,exact=True)*comb(N-6*i,x+j -6*i,exact=True)
    return val

def minproba(N,X,f): # f denotes the maximum number of faulty blocks allowed in the witness
    proba = 0 
    u = 0 
    val = [0 for k in range(N/6-1)] 
    for i in range(1,N/6+1): 
        for j in range(i+1): 
            u = Func_u(N,1,X,i,j) 
            val[0] += (-1)** (i+1)*u 
            for H in range(1,N/6-1): 
                if H != j-1: 
                    u = u*(H+1)*(N/6-H-i+j)/((N/6-H)*(H+1-j))
                else: 
                    u = Func_u(N,H+1,X,i,j)
                val[H] += (-1)** (i+1)*u
    
    for H in range(N/6-1): 
        val[H] = math.exp(math.log(val[H]) - math.log(comb(N,X, exact=True)))
    
    return min(val[f-1:N/6-f]) 

if __name__ == "__main__ ":
    f = 12 # number of faulty blocks allowed    
    K = 1776 # set the value of K here 
    d = 710 # indicates the number of values of theta to try, from K-d to K-( f+1).
           
        # to validate a parameter set , all values from 0.6*K to K-(f+1) must be verified
        # we conjecture that the best value is always K-(f+1)

for x in range(K-d, K-f): 
    p = 1-minproba(K,x,f) 
    if p==0.0: #if the failure probability of the adversary , given by minproba , is too close to 1 and Python rounds it to 1.0
        print("pmin("+str(x)+") ~= 0")
    else: 
        val = int(1/p) 
        print("pmin("+str(x)+") = 1/"+str(val))
