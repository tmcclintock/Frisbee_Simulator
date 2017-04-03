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

#The initial conditions
x,y,z = 0.0, 0.0, 1.0
vx,vy,vz = 10.0,0.0,0.0
phi,theta,gamma = 0.0,0.0,0.0
phidot,thetadot,gammadot = 0.0,0.0,50.0

#The frisbee
test_frisbee = frisbee.Frisbee(x,y,z,
                               vx,vy,vz,
                               phi,theta,gamma,
                               phidot,thetadot,gammadot)
model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])
test_frisbee.initialize_model(model)
coordinates = np.array([x,y,z,vx,vy,vz,phi,theta,gamma,\
                        phidot,thetadot,gammadot])

#Integrate it
trajectory = odeint(test_frisbee.equations_of_motion,coordinates,times)
x,y,z = trajectory[:,:3].T
err = np.ones_like(x)*0.05 #5 centimeters

#Plot it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
plt.plot(x,y,z)
plt.show()
plt.clf()

#Save it
outputs = np.array([times,x,y,z,err,err,err]).T
outputs = outputs[0::10] #Take only every ten
outputs = outputs[outputs[:,3]>0,:] #Take only entries with +z
np.savetxt("simulated_trajectory.txt",outputs,header="time (sec); x,y,z (m); x_err,y_err,z_err (m)")
