import pandas as pd
import numpy as np
import distances as dt
import channelModel as ch
import combinationsVectors as combV
import unimodularMatrix as benchmark

# System settings
Nt = 4  # Number of antennas at the transmitter
Pmax = 100  # Total available power
nUsers = 8  # Number of users
N = 2  # Number of subcarriers
dInnerRadius = 1
dOuterRadius = 10
gamma = 3  # Path loss exponent
num_iterations = 1000  # Number of repetitions
epsilon = 1e-4
uj = np.ones((nUsers, 1))  # Weight vector
Pn = [(Pmax / N) * x for x in np.ones((N, 1))]  # Initial power per subcarrier
# Channel generation
h = ch.channel(Nt, N, nUsers, gamma, dOuterRadius, dInnerRadius)
# User allocation algorithms
x_tum = benchmark.unimodularMatrixUserMatching(h, Pmax, N, nUsers, uj)
x_tum = x_tum.reshape((-1, 1))
# Vector combination
combVect = combV.combVector(nUsers, N)
userMatch = np.hstack((combVect, x_tum))

