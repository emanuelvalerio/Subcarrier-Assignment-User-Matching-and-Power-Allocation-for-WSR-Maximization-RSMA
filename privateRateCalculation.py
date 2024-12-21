import numpy as np
import normalizeVector as norm
def privateRate(h,n,ii,jj,Pn1,Pn2):
    hni = h[:,n,ii];
    hnj = h[:,n,jj]; 
    h_hni  =norm.normalization(h,n,ii);
    h_hnj  = norm.normalization(h,n,jj);
    rho = (1-(np.abs(np.dot(np.conj(h_hni).T , h_hnj))**2));
    gamma1 = 1+((np.linalg.norm(hni)**2)*rho*Pn1);
    gamma2 = 1+((np.linalg.norm(hnj)**2)*rho*Pn2);
    pRate1 = np.log2(gamma1);
    pRate2 = np.log2(gamma2);
    return pRate1,pRate2,gamma1,gamma2,rho