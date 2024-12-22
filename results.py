import pandas as pd
import matplotlib.pyplot as plt
import ast
# Passo 1: Carregar os dados
data = pd.read_csv("Results/monte_carlo_rates.csv")

# Passo 2: Converter strings para valores numéricos
# Aplica `ast.literal_eval` para transformar a string em lista e extrai o primeiro valor
for column in ["Rate_TUM", "Rate_Random", "Rate_EPA"]:
    data[column] = data[column].apply(lambda x: ast.literal_eval(x)[0] if isinstance(x, str) else x)

# Passo 3: Agrupar e calcular médias
mean_rates = data.groupby("N").mean().reset_index()

# Passo 4: Plotar o gráfico
plt.figure(figsize=(10, 6))
plt.plot(mean_rates["N"], mean_rates["Rate_TUM"], label="Rate_TUM", marker="o")
plt.plot(mean_rates["N"], mean_rates["Rate_Random"], label="Rate_Random", marker="s")
plt.plot(mean_rates["N"], mean_rates["Rate_EPA"], label="Rate_EPA", marker="^")

# Adicionar rótulos e título
plt.xlabel("Number of Subcarriers (N)")
plt.ylabel("Average Rate")
plt.title("Average Rate vs Number of Subcarriers")
plt.legend()
plt.grid(True)

# Exibir o gráfico
plt.show()
