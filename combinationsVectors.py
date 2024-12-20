import numpy as np
import math 
def combVector(nUsers,N):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    combV= np.zeros((N*comb,3)); 
    aux = 0;
    for n in range(0,N):
        prev = 0;
        for ii in range(0,comb):
            for jj in range( prev+1,nUsers):
                combV[aux,0] = int(n+1);
                combV[aux,1] = int(prev+1);
                combV[aux,2] = int(jj+1);
                aux += 1;
            prev += 1;
         
            if aux > comb*N:
                break;
    return combV;