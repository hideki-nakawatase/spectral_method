import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("K_dV_equation/data/fftw_linear_deleted_data.csv", header=None, sep=",")

x = df.iloc[:, 0]
y = df.iloc[:, 1]
y_sin = np.sin(x)

plt.plot(x, y)
plt.plot(x, y_sin)
plt.savefig("K_dV_equation/data/linear_delete.jpg")
