import numpy as np
import privateRateCalculation as prc
import commomRateCalculation as crc
import math
def partialDerivative():
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
    
def gradDes(Pmax,N):
    alpha = 15.1;
    Pn = [(Pmax / N) * x for x in np.ones((N, 1))]  # Initial power per subcarrier
    Pn_prev = Pn;
    t = 0;
    k = 0;
    z = [];
    
    while(t<1 or (np.abs(objectiveFun()-objectiveFun())>epsilon)):
        Pn_prev = Pn; # update the previous power vector
        dgk_dPk = partialDerivative();