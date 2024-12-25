import numpy as np
import math
import pandas as pd
import random

def randomUserMatching(combVect,N,nUsers):
    combVect_aux = np.hstack((combVect, np.zeros((combVect.shape[0],1))));
    ii = random.randint(1, nUsers);
    jj = random.randint(1, nUsers);
    while ii==jj:
       jj = random.randint(1, nUsers);
    
    max = jj;
    if ii>max:
        jj = ii;
        ii = max;
        
    prev = [];
    for n in range(N):
        nn = random.randint(1, N);
        while np.isin(nn,prev):
            nn = random.randint(1, N);
        condition = (combVect_aux[:, 0] == nn) & (combVect_aux[:, 1] == ii) & (combVect_aux[:, 2] == jj);
        combVect_aux[condition, 3] = 1
        prev.append(nn);
        
    return combVect_aux