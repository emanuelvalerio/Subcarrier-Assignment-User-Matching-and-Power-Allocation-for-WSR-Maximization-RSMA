import numpy as np
import pandas as pd
import distances as dt
def channel(Nt,N,nUsers,gamma,dOuterRadius,dInnerRadius):
   d = dt.distance(nUsers,dOuterRadius,dInnerRadius);
   h = np.zeros((Nt,N,nUsers), dtype=complex);
   for jj in range(0,nUsers):
      variance = 1/(Nt * d[jj]**gamma);
      for n in range(0,N):
         aux = (np.sqrt(variance/2))*(np.random.randn(1,Nt) + 1j*np.random.randn(1,Nt));
         h[:,n,jj] = aux;
   return h