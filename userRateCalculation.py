import numpy as np
import calculateFc
import normalizeVector as norm
import privateRateCalculation as prc
import commomRateCalculation as crc
import math

def userRate(h,restX,P_opt,uj,N,nUsers):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    userRate = np.zeros((nUsers,1));
    for n in range(N):
        aux = np.where(restX[:, 0].astype(int) == int(n+1))[0];
        idx = int(np.where(restX[aux[0]:aux[comb-1]+1,3] != 0)[0]);
        ii = int(restX[aux[idx],1]-1); # decrease in 1 because was increased 1 when the matrix was it built
        jj = int(restX[aux[idx],2]-1);
        hni = h[:,n,ii];
        hnj = h[:,n,jj];
        strongest = max(np.linalg.norm(hni),np.linalg.norm(hnj));
        if strongest == np.linalg.norm(hni):
            h1 = hni;
            h2 = hnj;
            u1 = uj[ii];
            u2 = uj[jj];
            h_1 = norm.normalization(h,n,ii);
            h_2 = norm.normalization(h,n,jj);
            rho = (1-(np.abs(np.dot(np.conj(h_1).T , h_2))**2));
            Pn1 = (1/3)*P_opt[n];
            Pn2 = (1/3)*P_opt[n];
            Pnc = (1/3)*P_opt[n];
            fc = calculateFc.calculationFc(h1,h2,Pn1,Pn2,rho);
            gamma1 = 1+(np.linalg.norm(h1)**2)*rho*Pn1;
            gamma2 = 1+(np.linalg.norm(h2)**2)*rho*Pn2;
        
            if(P_opt[n]!=0):
               userRate[ii] = userRate[ii] + math.log2(gamma1) + (1/2)*min(math.log2(1+(((np.abs(np.dot(np.conj(h2).T,fc))**2)*Pnc)/(gamma2))),np.log2(1+(((np.abs(np.dot(np.conj(h1).T,fc))**2)*Pnc)/(gamma1))));
               userRate[jj] = userRate[jj] + math.log2(gamma2) + (1/2)*min(math.log2(1+(((np.abs(np.dot(np.conj(h2).T,fc))**2)*Pnc)/(gamma2))),np.log2(1+(((np.abs(np.dot(np.conj(h1).T,fc))**2)*Pnc)/(gamma1))));
        else:
            h1 = hnj;
            h2 = hni;
            u1 = uj[jj];
            u2 = uj[ii];
            h_1 = norm.normalization(h,n,jj);
            h_2 = norm.normalization(h,n,ii);
            rho = (1-(np.abs(np.dot(np.conj(h_1).T , h_2))**2));
            Pn1 = (1/3)*P_opt[n];
            Pn2 = (1/3)*P_opt[n];
            Pnc = (1/3)*P_opt[n];
            fc = calculateFc.calculationFc(h1,h2,Pn1,Pn2,rho);
            gamma1 = 1+(np.linalg.norm(h1)**2)*rho*Pn1;
            gamma2 = 1+(np.linalg.norm(h2)**2)*rho*Pn2;
           
            if(P_opt[n]!=0):
               userRate[jj] = userRate[jj] + math.log2(gamma1) + (1/2)*min(math.log2(1+(((np.abs(np.dot(np.conj(h2).T,fc))**2)*Pnc)/(gamma2))),np.log2(1+(((np.abs(np.dot(np.conj(h1).T,fc))**2)*Pnc)/(gamma1))));
               userRate[ii] = userRate[ii] + math.log2(gamma2) + (1/2)*min(math.log2(1+(((np.abs(np.dot(np.conj(h2).T,fc))**2)*Pnc)/(gamma2))),np.log2(1+(((np.abs(np.dot(np.conj(h1).T,fc))**2)*Pnc)/(gamma1))));
    return userRate