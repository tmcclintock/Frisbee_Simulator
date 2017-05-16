import numpy as np
import corner
import matplotlib.pyplot as plt


chainpath = "chains/PD0PDa_chain.txt"
fullchain = np.loadtxt(chainpath)
nburn = len(fullchain)/10 #take 10%
chain = fullchain[nburn:]
fig = corner.corner(chain)
plt.show()
