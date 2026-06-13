import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

df = pd.read_csv(
    "K_dV_equation/data/fftw_linear_deleted_data_time_step.csv",
    header=None,
    sep=",",
)
df1 = pd.read_csv("K_dV_equation/data/fftw_data.csv", header=None, sep=",")

fig, ax = plt.subplots()
ims = []
x = df1.iloc[:, 0]

ax.set_xlim(x.min(), x.max())
ax.set_ylim(-2, 6)
for i in range(len(df)):
    y = df.iloc[i, :]
    ax.set_title(f"step={i*100}")
    im = ax.plot(x, y, color="blue")
    ims.append(im)

ani = animation.ArtistAnimation(fig, ims, interval=100)
ani.save("K_dV_equation/data/linear_delete.gif", writer="pillow")
plt.close()
print("animation done")
