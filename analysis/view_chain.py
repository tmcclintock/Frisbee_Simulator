import numpy as np
import corner
import matplotlib.pyplot as plt

nburn = int(4000)
fullchain = np.loadtxt("Final_PL0,PLa,PD0,PDa_nsteps=20,000_nwalkers=8.txt")
chain = fullchain[nburn:]
fig = corner.corner(chain)
#plt.xlim(0.32,0.34)
plt.show()
