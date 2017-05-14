"""
Create figure 1, which shows the trajectory. Create this both with and without scatter.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True, fontsize=24)
plt.rc("errorbar", capsize=3)

datapath = "../simulation_data/sample_throw.txt"
t, x, y, z, xe, ye, ze = np.loadtxt(datapath, unpack=True)
xs = x + np.random.randn(len(x))*0.05
ys = y + np.random.randn(len(y))*0.05
zs = z + np.random.randn(len(z))*0.05

#original throw was integrated at 300 fps, or dt = 1/300.
#downsampling by 10 mean we have a 30 fps throw, which is realistic for most cameras
ds = 10
plt.plot(x, z)
plt.plot(xs[::ds], zs[::ds], zorder=-1)
plt.xlabel(r"$x$ displacement $[{\rm m}]$")
plt.ylabel(r"$z$ displacement $[{\rm m}]$")
plt.subplots_adjust(bottom=0.17, left=0.15)
plt.gcf().savefig("figures/xz_scatter.pdf")
plt.show()
