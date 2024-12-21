import numpy as np
import math
import random

def randomUserMatching(N,nUsers):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    x = [];
    for n in range(N):
        aux = np.zeros((comb,1));
        idx_sorted = random.randint(1, comb);
        aux[idx_sorted] = 1;
        x.append(aux);
    return x;