"""
Compare the output of getting the trajectory from the C code with the python code.
"""
import sys
sys.path.insert(0,"../")
import frisbee
import numpy as np

#Set up an array of times
time_initial = 0.0
time_final = 3.0
frac = 1000
dt = 0.09999/frac
N_times = int((time_final-time_initial)/dt) #300 FPS
times = np.linspace(time_initial,time_final,N_times)

#The initial conditions
x,y,z = 0.0, 0.0, 1.0
vx,vy,vz = 10.0,0.1,0.1
phi,theta,gamma = 0.1,0.1,0.1
phidot,thetadot,gammadot = 0.1,0.1,50.0

#The python frisbee
py_frisbee = frisbee.Frisbee(x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot)
model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])
py_frisbee.initialize_model(model)
ptimes,ptraj = py_frisbee.get_trajectory(time_initial, time_final, dt)

#The c frisbee
c_frisbee = frisbee.Frisbee(x,y,z,vx,vy,vz,phi,theta,gamma,phidot,thetadot,gammadot,use_C=True)
model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])
c_frisbee.initialize_model(model)
ctimes,ctraj = c_frisbee.get_trajectory(time_initial, time_final, dt)

xp = ptraj[:,0]
xc = ctraj[:,0]
yp = ptraj[:,1]
yc = ctraj[:,1]
zp = ptraj[:,2]
zc = ctraj[:,2]
phip = ptraj[:,6]
phic = ctraj[:,6]

for i in range(len(phip)):
    print zp[i*frac], zc[i*frac]
