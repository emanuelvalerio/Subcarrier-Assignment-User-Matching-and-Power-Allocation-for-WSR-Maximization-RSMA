import numpy as np
import calculateFc
import normalizeVector as norm
import privateRateCalculation as prc
import commomRateCalculation as crc
import math

def userRate(h,restX,P_opt,uj,N,nUsers):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    rate_subcarrier = np.zeros((N,1));
    user_rate=np.zeros((nUsers,1));
    for n in range(N):
        aux = np.where(restX[:, 0].astype(int) == int(n+1))[0];
        idx = int(np.where(restX[aux[0]:aux[comb-1]+1,3] != 0)[0]);
        ii = int(restX[aux[idx],1]-1); # decrease in 1 because was increased 1 when the matrix was it built
        jj = int(restX[aux[idx],2]-1);
        hni = h[:,n,ii];
        hnj = h[:,n,jj];
        h_ni = norm.normalization(h,n,ii);
        h_nj = norm.normalization(h,n,jj);
        rho = (1-(np.abs(np.dot(np.conj(h_ni).T , h_nj))**2));
        Pn1 = (1/3)*P_opt[n];
        Pn2 = (1/3)*P_opt[n];
        Pnc = (1/3)*P_opt[n];
        fc = calculateFc.calculationFc(hni,hnj,Pn1,Pn2,rho);
        gamma1 = 1+(np.linalg.norm(hni)**2)*rho*Pn1;
        gamma2 = 1+(np.linalg.norm(hnj)**2)*rho*Pn2;
        pRate1 = np.log2(gamma1);
        pRate2 = np.log2(gamma2);
        cRate = np.log2(1+(((np.abs(np.dot(np.conj(hnj).T,fc))**2)*Pnc)/gamma2));
        rate_subcarrier[n] = pRate1+pRate2+cRate;
        user_rate[ii] = user_rate[ii]+pRate1+(1/2)*cRate;
        user_rate[jj] = user_rate[jj]+pRate2+(1/2)*cRate;
    return rate_subcarrier,user_rate;