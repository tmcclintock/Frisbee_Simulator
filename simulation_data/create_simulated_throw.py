"""
This is a script that creates a simulated throw that we 
will then analyse.
"""
import sys
sys.path.insert(0,"../")
import frisbee
import numpy as np
from scipy.integrate import odeint

#Set up an array of times
time_initial = 0.0
time_final = 3.0
N_times = int((time_final-time_initial)/0.003333) #300 FPS
times = np.linspace(time_initial,time_final,N_times)
dt = times[1]-times[0]

#The initial conditions
x,y,z = 0.0, 0.0, 1.0
vx,vy,vz = 10.0,0.0,0.0
phi,theta,gamma = 0.0,-0.25,0.0
phidot,thetadot,gammadot = 0.0,0.0,50.0

#The frisbee
test_frisbee = frisbee.Frisbee(x,y,z,
                               vx,vy,vz,
                               phi,theta,gamma,
                               phidot,thetadot,gammadot,
                               use_C=True)
model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])
test_frisbee.initialize_model(model)
times, trajectory = test_frisbee.get_trajectory(time_initial, time_final, dt=dt)
print trajectory.shape
x, y, z = trajectory.T[:3]

#Create another throw with 5 cm scatter in each direction
xs = x + 0.05*np.random.randn(len(x))
ys = y + 0.05*np.random.randn(len(x))
zs = z + 0.05*np.random.randn(len(x))

#Plot it and the scattered throw
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
plt.plot(x, y, z)
plt.plot(xs, ys, zs, alpha=0.5)
plt.show()
plt.clf()

#Save it
err = np.ones_like(x) * 0.05 #5 centimeter error
outputs = np.array([times, x, y, z, err, err, err]).T
outputs = outputs[0::1] #Take every one
outputs = outputs[outputs[:,3]>0,:] #Take only entries with +z
np.savetxt("sample_throw.txt",outputs,header="time (sec); x,y,z (m); x_err,y_err,z_err (m)")

#Now save the full trajectory
fulloutputs = np.vstack((times,trajectory.T)).T
fulloutputs = fulloutputs[fulloutputs[:,3]>0,:] #Take only entries with +z
np.savetxt("simulated_trajectory.txt",fulloutputs,header="t x y z vx vy vz phi theta gamma phid thetad gammad")
