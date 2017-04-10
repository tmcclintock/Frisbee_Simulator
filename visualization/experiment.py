"""
This is an in-development vesion of plot_trajectory.py.

The plan:
show an animation of the trajectory
also have the pitch and rolls in adjacent plots
also have the spin vs time in an adjacent plot
"""

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import gridspec

def update_lines(num, dataLines, lines):
    """
    #template for lines:
    for i in range(len(lines)):
    line = lines[i]
    data = dataLines[i]
    switch(i):
    case 0: 3D line
    case 1: polar plot
    case 2: 2D line
    return lines
    """

    #Each dataLine is a 3xnum array
    #lines contains the matplotlib line objects
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2, :num])
        if len(data) == 3:
            line.set_3d_properties(data[2, :num])
    return lines

def update_line(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
    return line

trajectory=np.loadtxt("full_trajectory.txt")
print trajectory.shape
trajectory = trajectory[::4] #Reduce the size for the gif purposes

times,x,y,z = trajectory.T[:4]
print trajectory.shape
phi, theta, gamma = trajectory.T[7:10]
phid, thetad, gammad = trajectory.T[10:]
print gamma[0]

zeros = np.zeros(len(x))

fig = plt.figure()
#The gridspec shit
gs = gridspec.GridSpec(3, 3)
gs.update(left=0.05, right=0.78, wspace=0.55)

#These are the axes and lines for the 3D plot
ax0 = plt.subplot(gs[:, :2], projection='3d')
ax0.set_title('Trajectory')
ax0.set_zlim(0,max(z))
line = ax0.plot(x, y ,z)[0]
wall_shadow   = ax0.plot(x, zeros+max(y), z, linestyle='--')[0]
ground_shadow = ax0.plot(x, y, zeros, linestyle='--')[0]

#While the 3D plot contains the trajectory,
#these extra plots can contain the pitch and roll angles, updated
#along with the trajectory
ax1 = plt.subplot(gs[0,2], projection='polar')
ax1.set_rmax(1)
ax1.set_rticks([])
philine = ax1.plot(phi,times/max(times))[0]
ax2 = plt.subplot(gs[1,2], projection='polar')
ax2.set_rmax(1)
ax2.set_rticks([])
thetaline = ax2.plot(theta,times/max(times))[0]
ax3 = plt.subplot(gs[2,2])
gammadline = ax3.plot(times, gammad)[0]


dataLines = [[x, y, z], [x, zeros+max(y), z], [x, y, zeros], [phi, times/max(times)], [theta, times/max(times)], [times, gammad]]
dataLines = [np.array(DL) for DL in dataLines] #Gotta convert to numpy arrays
lines = [line, wall_shadow, ground_shadow, philine, thetaline, gammadline]
print len(x)," frames in the animation"
print "This will take %f milliseconds"%(len(x)*0.5)

anim = animation.FuncAnimation(fig, update_lines, frames=len(x), 
                               fargs=(dataLines, lines), interval=1, blit=True)

plt.show()
plt.clf()
