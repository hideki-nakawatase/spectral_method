import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("non_linear/spectral_non_linear.csv", header=None, sep=",")

x = df.iloc[:, 0]
y1 = df.iloc[:, 1]

plt.figure(figsize=(12, 8))
plt.plot(x, y1)
plt.show()
