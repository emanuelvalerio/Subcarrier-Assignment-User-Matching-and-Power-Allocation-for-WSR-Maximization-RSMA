import pandas as pd
import numpy as np
import distances as dt
import channelModel as ch
import lowComplexityUserMatching as uMatch
import combinationsVectors as combV
import optimalPA as optPA
import rateCalculation
import unimodularMatrixBenchmark as benchmark

# Configurações do sistema
Nt = 4  # Número de antenas no transmissor
Pmax = 100  # Potência total disponível
nUsers = 8 # Número de usuários
dInnerRadius = 1
dOuterRadius = 10
gamma = 3  # Expoente de perda de caminho
num_iterations = 1000  # Número de repetições
epsilon = 1e-4
uj = np.ones((nUsers, 1))  # Vetor de pesos
Pn = [(Pmax / N) * x for x in np.ones((N, 1))]  # Potência inicial por subportadora
# Geração do canal
h = ch.channel(Nt, N, nUsers, gamma, dOuterRadius, dInnerRadius)
 # Algoritmos de alocação de usuários
x_tum = benchmark.unimodularMatrixUserMatching(h, Pmax, N, nUsers, uj)
x_tum = x_tum.reshape((-1,1));
# Combinação de vetores
combVect = combV.combVector(nUsers, N)
userMatch = np.hstack((combVect, x_tum))
