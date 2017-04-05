import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,

trajectory=np.loadtxt("sample_throw.txt")

x = trajectory[:,0]
y = trajectory[:,1]
z = trajectory[:,2]

zeros = np.zeros(len(x))

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

line = ax.plot(x, y ,z)
#shadow_1 = ax.plot(zeros, y, z)
#shadow_2 = ax.plot(x, zeros, z)

#plt.plot(x,y,z)
plt.show()
plt.clf()
