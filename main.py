import pandas as pd
import numpy as np
import distances as dt
import channelModel as ch
import combinationsVectors as combV
import unimodularMatrix as tum
import gradientDescent as gradDes
import randomMatching as randMatch
import userRateCalculation

# System settings
Nt = 4  # Number of antennas at the transmitter
Pmax = 100  # Total available power
nUsers = 8  # Number of users
dInnerRadius = 1
dOuterRadius = 10
gamma = 3  # Path loss exponent
num_iterations = 1000  # Number of repetitions
epsilon = 1e-4
uj = np.ones((nUsers, 1))  # Weight vector
n_iterations = 1000  # Número de repetições de Monte Carlo

N_values = [10];  # Number of subcarrier that will be evaluation

results = []

for N in N_values:
    print(N)
    Pn = [(Pmax / N) * x for x in np.ones((N, 1))]  # Initial power per subcarrier
    for i in range(n_iterations):
        print("Iteration " + str(i))

        # Channel Generation
        h = ch.channel(Nt, N, nUsers, gamma, dOuterRadius, dInnerRadius)
        
        # User Matching algorithms
        combVect = combV.combVector(nUsers, N)
        x_tum = tum.unimodularMatrixUserMatching(h, Pmax, N, nUsers, uj)
        x_tum = x_tum.reshape((-1, 1))
        userMatch_rand = randMatch.randomUserMatching(combVect,N, nUsers)
        
        userMatch_tum = np.hstack((combVect, x_tum))
        
        # Power Optimization
        Pn_opt_tum = gradDes.gradDes(h, userMatch_tum, nUsers, Pmax, N, uj, epsilon);
        Pn_opt_rand = gradDes.gradDes(h, userMatch_rand, nUsers, Pmax, N, uj, epsilon)
        
        # Cálculo da taxa total
        rate_tum,user_rate_tum= userRateCalculation.userRate(h,userMatch_tum,Pn_opt_tum,uj,N,nUsers);
        rate_rand,user_rate_rand= userRateCalculation.userRate(h,userMatch_rand,Pn_opt_rand,uj,N,nUsers);
        rate_EPA,user_rate_EPA=userRateCalculation.userRate(h,userMatch_tum,Pn,uj,N,nUsers);
        
        # Salvar os resultados na lista
        results.append([N, i + 1, sum(user_rate_tum).item(), sum(user_rate_rand).item(),sum(user_rate_EPA).item()])

# Converter os resultados em um DataFrame do pandas
df = pd.DataFrame(results, columns=["N", "Repetition", "Rate_TUM", "Rate_Random","Rate_EPA"])

# Salvar os resultados em um arquivo CSV
df.to_csv("Results/monte_carlo_rates_aux.csv", index=False)

print("Simulação Monte Carlo finalizada. Resultados salvos em 'monte_carlo_rates.csv'.")
