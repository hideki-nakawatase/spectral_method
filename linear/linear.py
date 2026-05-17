import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("linear/linear.csv")

N = 100
x = df.iloc[:,0]
y1 = df.iloc[:,1]
y2 = df.iloc[:,2]

plt.figure(figsize=(12,8))
plt.plot(x,y1)
plt.plot(x,y2,ls=":")
plt.show()