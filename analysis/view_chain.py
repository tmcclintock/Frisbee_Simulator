import numpy as np
import corner
import matplotlib.pyplot as plt

nburn = int(100)
fullchain = np.loadtxt("fullchain.txt")
chain = fullchain[nburn:]
fig = corner.corner(chain)
plt.xlim(0.32,0.34)
plt.show()
