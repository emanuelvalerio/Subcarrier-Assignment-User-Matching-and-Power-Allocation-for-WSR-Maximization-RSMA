import pandas as pd
import matplotlib.pyplot as plt
import ast
# Passo 1: Carregar os dados
data = pd.read_csv("Results/monte_carlo_rates_aux.csv")

# Passo 2: Converter strings para valores numéricos
# Aplica `ast.literal_eval` para transformar a string em lista e extrai o primeiro valor
for column in ["Rate_TUM", "Rate_Random", "Rate_EPA"]:
    data[column] = data[column].apply(lambda x: ast.literal_eval(x)[0] if isinstance(x, str) else x)

# Passo 3: Agrupar e calcular médias
mean_rates = data.groupby("N").mean().reset_index()

# Passo 4: Plotar o gráfico
plt.figure()
plt.plot(mean_rates["N"], mean_rates["Rate_TUM"], label="Gradient Descent", marker="^",color='r')
plt.plot(mean_rates["N"], mean_rates["Rate_EPA"], label="EPA", marker="o",color='b')
plt.plot(mean_rates["N"], mean_rates["Rate_Random"], label="Random Matching", marker="s",color='g')


# Adicionar rótulos e título
plt.xlabel("Number of Subcarriers (N)")
plt.ylabel("Sum-rate (bps/Hz)")
plt.title("Sum-Rate vs Number of Subcarriers")
plt.legend(loc="lower right")  # Posição predefinida
plt.grid(True, linestyle='-', linewidth=0.3, color='gray')
plt.xticks([4,5,6,7,8,9,10])
plt.yticks([0,0.4,0.8,1.2,1.6,2.0,2.4,2.8])
# Exibir o gráfico
plt.show()
