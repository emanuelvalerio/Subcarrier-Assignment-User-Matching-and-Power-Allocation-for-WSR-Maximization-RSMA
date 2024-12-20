import numpy as np

def normalization(h,n,idx):
    return h[:,n,idx]/np.linalg.norm(h[:,n,idx]);
    