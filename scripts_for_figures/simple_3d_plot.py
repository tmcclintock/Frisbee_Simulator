import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def update_lines(num, dataLines, lines):
    #Each dataLine is a 3xnum array
    #lines contains the matplotlib line objects
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

def update_line(num, data, line):
    line.set_data(data[:2, :num])
    line.set_3d_properties(data[2, :num])
    return line

#Read in the data and down-sample to get the FPS we want
trajectory = np.loadtxt("../simulation_data/sample_throw.txt")
T = trajectory[-1, 0] -trajectory[0, 0] #Total time
N = len(trajectory)
FPS = 60
trajectory = trajectory[::int(N/(FPS*T))]

times,x,y,z = trajectory.T[:4]

zeros = np.zeros(len(x))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Create the line for the trajectory
line = ax.plot(x, y ,z)[0]

#Create a shadow on the far (y_max) wall and on the ground
wall_shadow   = ax.plot(x, zeros+max(y), z, linestyle='--')[0]
ground_shadow = ax.plot(x, y, zeros, linestyle='--')[0]
dataLines = [[x, y, z], [x, zeros+max(y), z], [x, y, zeros]]
dataLines = [np.array(DL) for DL in dataLines] #Gotta convert to numpy arrays
lines = [line,wall_shadow,ground_shadow]

print len(x)," frames in the animation"
print "This will take %f milliseconds"%(len(x)*0.5)
#anim = animation.FuncAnimation(fig, update_lines, frames=len(x), fargs=(dataLines, lines), interval=5, blit=True)
#anim.save('line.gif', dpi=80, writer='imagemagick', fps=60)

plt.show()
plt.clf()
