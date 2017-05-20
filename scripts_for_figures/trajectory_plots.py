"""
Create figures showing the trajectory, and trajectories with scatter on them.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rc("text", usetex=True, fontsize=24)
plt.rc("errorbar", capsize=3)

def traj_2D(x,y,z,xs,ys,zs):
    plt.plot(x, z, c='b')
    plt.plot(xs, zs, zorder=-1, c='k')
    plt.xlabel(r"$x\ {\rm displacement}\ [{\rm m}]$")
    plt.ylabel(r"$z\ {\rm displacement}\ [{\rm m}]$")
    plt.subplots_adjust(bottom=0.17, left=0.15)
    plt.gcf().savefig("figures/xz_scatter.pdf")
    plt.show()

def traj_2D_with_err(x,y,z,xs,ys,zs):
    plt.plot(x, z, c='b')
    plt.errorbar(xs, zs, xerr=xe, yerr=ze, zorder=-1, c='k', ls='')
    plt.xlabel(r"$x\ {\rm displacement}\ [{\rm m}]$")
    plt.ylabel(r"$z\ {\rm displacement}\ [{\rm m}]$")
    plt.subplots_adjust(bottom=0.17, left=0.15)
    plt.gcf().savefig("figures/xz_scatter_error.pdf")
    plt.show()

def traj_3D(x,y,z,xs,ys,zs):
    plt.rc("text", usetex=True, fontsize=16)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    plt.plot(x, y, z, c='b')
    plt.plot(xs, ys, zs, alpha=0.5, c='k')
    plt.xlabel(r"$x\ {\rm displacement}\ [{\rm m}]$")
    plt.ylabel(r"$y\ {\rm displacement}\ [{\rm m}]$")
    ax.set_zlabel(r"$z\ {\rm displacement}\ [{\rm m}]$")
    plt.gcf().savefig("figures/3d_scatter.pdf")
    plt.show()

if __name__ == "__main__":
    #Read in the data and down-sample to get the FPS we want
    datapath = "../simulation_data/sample_throw.txt"
    trajectory = np.loadtxt(datapath)
    T = trajectory[-1, 0] -trajectory[0, 0] #Total time
    N = len(trajectory)
    FPS = 60
    trajectory = trajectory[::int(N/(FPS*T))]
    print trajectory.shape
    t, x, y, z, xe, ye, ze = trajectory.T
    xs = x + np.random.randn(len(x))*0.05
    ys = y + np.random.randn(len(y))*0.05
    zs = z + np.random.randn(len(z))*0.05

    traj_2D(x,y,z,xs,ys,zs)
    traj_2D_with_err(x,y,z,xs,ys,zs)
    traj_3D(x,y,z,xs,ys,zs)
