import numpy as np
import privateRateCalculation as prc
import commomRateCalculation as crc
import math

def partialDerivative(h,restX,Pn,nUsers,n,uj):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    aux = np.where(restX[:, 0].astype(int) == int(n+1))[0];
    idx = int(np.where(restX[aux[0]:aux[comb-1]+1,3] != 0)[0]);
    ii = int(restX[aux[idx],1]-1); # decrease in 1 because was increased 1 when the matrix was it built
    jj = int(restX[aux[idx],2]-1);
    Pn1 = Pn/3;
    Pn2 = Pn/3;
    Pnc = Pn/3;
    hni = h[:,n,ii];
    hnj = h[:,n,jj];
    pRate1,pRate2,gamma1,gamma2,rho = prc.privateRate(h,n,ii,jj,Pn1,Pn2);
    Nnij,Dnij,cRate = crc.commomRate(h,Pnc,n,ii,jj,gamma1,gamma2,rho);
    dRnij1_dPn = (1/np.log(2))*(((np.linalg.norm(hni)**2)*rho/3)/(1+((np.linalg.norm(hni)**2)*rho*Pn1)));
    dRnij2_dPn = (1/np.log(2))*(((np.linalg.norm(hnj)**2)*rho/3)/(1+((np.linalg.norm(hnj)**2)*rho*Pn2)));
    auxDerivative = (np.linalg.norm(hni)**2) + (np.linalg.norm(hnj)**2) - (((gamma1+gamma2)/(np.sqrt(gamma1)*np.sqrt(gamma2)))*(np.abs(np.dot(np.conj(hni).T,hnj))**2));
    dRnijc_dPn = (1/np.log(2))*((rho*(np.linalg.norm(hni)**2)* (np.linalg.norm(hnj)**2)*auxDerivative)/(3*Dnij*(Nnij+Dnij)));
    df_dPn = uj[ii]*dRnij1_dPn + uj[jj]*dRnij2_dPn + ((uj[ii]+uj[jj])/2)*dRnijc_dPn;
    return df_dPn;

def objectiveFun(h,restX,Pn,N,nUsers,uj):
    comb = int((math.factorial(nUsers)/((math.factorial(nUsers-2))*math.factorial(2))));
    faux = np.zeros((N,1));
    
    for n in range(N):
        aux = np.where(restX[:, 0].astype(int) == int(n+1))[0];
        idx = int(np.where(restX[aux[0]:aux[comb-1]+1,3] != 0)[0]);
        ii = int(restX[aux[idx],1]-1); # decrease in 1 because was increased 1 when the matrix was it built
        jj = int(restX[aux[idx],2]-1);
        Pn1 = Pn2 = Pnc = Pn[n]/3;
        pRate1,pRate2,gamma1,gamma2,rho = prc.privateRate(h,n,ii,jj,Pn1,Pn2);
        Nnin,Dnij,cRate = crc.commomRate(h,Pnc,n,ii,jj,gamma1,gamma2,rho);
        faux[n,:] = uj[ii]*pRate1 + uj[jj]*pRate2 + ((uj[ii]+uj[jj])/2)*cRate;
        
    return sum(faux);

def gradDes(h,restX,nUsers,Pmax,N,uj,epsilon):
    alpha = 7.1; # step size
    Pn = [(Pmax / N) * x for x in np.ones((N, 1))]  # Initial power per subcarrier
    Pn_prev = [(Pmax / N) * x for x in np.ones((N, 1))];
    t = 0;
    k = 0;
    Z = [];
    
    while(t<1 or (np.abs(objectiveFun(h,restX,Pn,N,nUsers,uj)-objectiveFun(h,restX,Pn_prev,N,nUsers,uj))>epsilon)):
        Pn_prev = Pn.copy(); # update the previous power vector
        dgk_dPk = partialDerivative(h,restX,Pn[k],nUsers,k,uj);
        for n in set(range(N)) - set([k] + Z):
            #Update the power
            dgn_dPn = partialDerivative(h,restX,Pn[n],nUsers,n,uj);
            dfk_dPn = dgn_dPn - dgk_dPk;
            Pn[n] = Pn_prev[n] + alpha*dfk_dPn;
            
        for ii in range(len(Pn)):
            if Pn[ii]<0:
                Z.extend([ii]);# Find the indexs of negative power allocated
                
        for z in range(len(Z)):
            Pn[Z[z]] = 0; # Power constraint, so, for all negative power allocated we replaced to zero.
        
        if np.isin(k,Z):
            aux = list(set(range(0, N)) - set(Z));
            k = aux[0]; # update k
            
        #Upper Bound Power constraint
        Pn[k] = Pmax - sum(Pn[i] for i in set(range(N)) - {k});
        
        t+=1;
       # alpha=alpha/t;
        
    return Pn;
        
