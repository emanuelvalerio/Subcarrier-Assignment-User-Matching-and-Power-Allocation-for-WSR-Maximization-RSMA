import numpy as np
import math
import pandas as pd
import random

def randomUserMatching(combVect,N,nUsers):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    x = [];
    for n in range(N):
        auxMatrix = combVect[n*comb:(n+1)*comb,:].copy();
        np.random.shuffle(auxMatrix);
        aux = np.zeros((comb,1));
        idx_sorted = random.randint(0, comb-1);
        aux[idx_sorted] = 1;
        x.extend( np.hstack((auxMatrix, aux)));
    df = pd.DataFrame(x) # Convert the numpy array to a pandas DataFrame
    return df.sort_values(by=[0, 1, 2]).values # Sort by the first, second, and third columns
    
