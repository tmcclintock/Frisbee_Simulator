"""
Create the multi-panel gif with the angles displayed in polar plots.
"""

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib import gridspec
plt.rc('text',usetex=True)

def update_lines(num, dataLines, lines):
    for i in range(len(lines)):
        line = lines[i]
        data = dataLines[i]
        if i < 3: #Trajectory plot
            line.set_data(data[0:2, :num])
            line.set_3d_properties(data[2, :num])
        elif 2 < i <5: #Angle plots
            line.set_data([data[0][num],data[0][num]+np.pi],[1,1])
        elif i == 5: #Spin vs time
            line.set_data(data[0:2, :num])

    return lines
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

#Read in the data and down-sample to get the FPS we want
trajectory=np.loadtxt("../simulation_data/full_trajectory.txt")
T = trajectory[-1, 0] -trajectory[0, 0] #Total time
N = len(trajectory)
FPS = 60
trajectory = trajectory[::int(N/(FPS*T))]

times,x,y,z = trajectory.T[:4]
phi, theta, gamma = trajectory.T[7:10]
phi *= -1 #This is to orient things correctly for the plot
theta *= -1 #This is to orient things correctly for the plot
phid, thetad, gammad = trajectory.T[10:]

zeros = np.zeros(len(x))

fig = plt.figure()
gs = gridspec.GridSpec(3, 3)
gs.update(left=0.05, right=0.78, wspace=0.55)

#These are the axes and lines for the 3D plot
ax0 = plt.subplot(gs[:, :2], projection='3d')
ax0.set_title(r'${\rm Trajectory}$')
ax0.set_zlim(0,max(z))
line = ax0.plot(x, y ,z)[0]
wall_shadow   = ax0.plot(x, zeros+max(y), z, linestyle='--')[0]
ground_shadow = ax0.plot(x, y, zeros, linestyle='--')[0]
ax0.set_xlabel(r"${\rm X}\ [{\rm m}]$")
ax0.set_ylabel(r"${\rm Y}\ [{\rm m}]$")
ax0.set_zlabel(r"${\rm Z}\ [{\rm m}]$")

#While the 3D plot contains the trajectory,
#these extra plots can contain the pitch and roll angles, updated
#along with the trajectory
title_loc = 0.98
ax1 = plt.subplot(gs[0,2], projection='polar')
philine = ax1.plot([phi[0],phi[0]+np.pi], [1,1])[0]
ax1.set_rmax(1)
ax1.set_rticks([])
ax1.set_title(r"${\rm Left-Right}$",y=title_loc)
ax1.set_xticklabels([])

ax2 = plt.subplot(gs[1,2], projection='polar')
thetaline = ax2.plot([theta[0],theta[0]+np.pi],[1,1])[0]
ax2.set_rmax(1)
ax2.set_rticks([])
ax2.set_xticklabels([])
ax2.set_title(r"${\rm Front-Back}$",y=title_loc)

ax3 = plt.subplot(gs[2,2])
gammadline = ax3.plot(times, gammad)[0]
ax3.set_ylabel(r"${\rm Spin}$ $[{\rm rad/s}]$")
ax3.set_xlabel(r"${\rm Time}$ $[{\rm s}]$")
ax3.set_title(r"${\rm Rotation\ Rate}$",y=title_loc)


dataLines = [[x, y, z], [x, zeros+max(y), z], [x, y, zeros], [phi], [theta], [times, gammad]]
dataLines = [np.array(DL) for DL in dataLines] #Gotta convert to numpy arrays
lines = [line, wall_shadow, ground_shadow, philine, thetaline, gammadline]
print len(x)," frames in the animation"
print "This will take %f milliseconds"%(len(x)*0.5)

anim = animation.FuncAnimation(fig, update_lines, frames=len(x), 
                               fargs=(dataLines, lines), interval=10, blit=True)

#anim.save('multi_plot.gif', dpi=80, writer='imagemagick', fps=60)

plt.show()
plt.clf()
