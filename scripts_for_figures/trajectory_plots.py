"""
Create figure 1, which shows the trajectory. Create this both with and without scatter.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rc("text", usetex=True, fontsize=24)
plt.rc("errorbar", capsize=3)

datapath = "../simulation_data/sample_throw.txt"
t, x, y, z, xe, ye, ze = np.loadtxt(datapath, unpack=True)
xs = x + np.random.randn(len(x))*0.05
ys = y + np.random.randn(len(y))*0.05
zs = z + np.random.randn(len(z))*0.05



def traj_2D():
    ds = 10
    plt.plot(x, z, c='b')
    plt.plot(xs[::ds], zs[::ds], zorder=-1, c='k')
    plt.xlabel(r"$x\ {\rm displacement}\ [{\rm m}]$")
    plt.ylabel(r"$z\ {\rm displacement}\ [{\rm m}]$")
    plt.subplots_adjust(bottom=0.17, left=0.15)
    plt.gcf().savefig("figures/xz_scatter.pdf")
    plt.show()

def traj_3D():
    ds = 10
    plt.rc("text", usetex=True, fontsize=16)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    plt.plot(x, y, z, c='b')
    plt.plot(xs[::ds], ys[::ds], zs[::ds], alpha=0.5, c='k')
    plt.xlabel(r"$x\ {\rm displacement}\ [{\rm m}]$")
    plt.ylabel(r"$y\ {\rm displacement}\ [{\rm m}]$")
    ax.set_zlabel(r"$z\ {\rm displacement}\ [{\rm m}]$")
    plt.gcf().savefig("figures/3d_scatter.pdf")
    plt.show()

if __name__ == "__main__":
    traj_2D()
    traj_3D()
