import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("non_linear/conversion_method_fftw.csv", header=None, sep=",")

x = df.iloc[:, 0]
y1 = df.iloc[:, 1]

plt.figure(figsize=(12, 8))
plt.plot(x, y1)
plt.savefig("non_linear/conversion_method_fftw_t_1.05.jpg")
plt.show()
