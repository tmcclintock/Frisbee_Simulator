"""
Perform our analysis, but subsample the trajectory data.
This essentially recreates any FPS camera we choose.
"""
import sys
import numpy as np
import emcee, corner
import matplotlib.pyplot as plt
from analysis_routines import lnprob

#Read in the trajectory 
data = np.genfromtxt("../simulation_data/sample_throw.txt").T

#True starting condtions of kinematic variables we don't test
xi, yi, zi                 = data[1:4, 0] #meters
vxi, vyi, vzi              = 10.0, 0.0, 0.0 #m/s
phi, theta, gamma          = 0.0, -0.25, 0.0 #radians
phidot, thetadot, gammadot = 0.0, 0.0, 50.0 #rad/sec
initial_conditions = [xi, yi, zi, 
                      vxi, vyi, vzi, 
                      phi, theta, gamma, 
                      phidot, thetadot, gammadot]

#Figure out the initial FPS of the data
t = data[0, :]
N = len(t)
FPSi = int(round(1/(t[1]-t[0])))
print "Initial 'FPS':",FPSi

#Loop from FPSi to having only 3 frames through the whole trajectory
last_FPS = 999999
for i in range(1, N/3):
    ti = t[::i]
    FPScurrent = int(round(1/(ti[1]-ti[0])))
    if FPScurrent == last_FPS: 
        continue
    last_FPS = FPScurrent
    print "%d data points"%len(ti), "at FPS = %.2f"%(1./(ti[1]-ti[0]))
    #subsample the data
    datai = data[:, ::i]
    
    #Set up the walkers in parameter space with true positions
    test_params = [0.18, 0.69] #PD0, PDa true positions

    nwalkers = 4
    ndim = 2
    nsteps = 500
    pos = [test_params + 1e-2*np.fabs(test_params)*np.random.randn(ndim) 
           for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(datai, initial_conditions))
    sampler.run_mcmc(pos,nsteps)
    fullchain = sampler.flatchain
    chain = fullchain[200:]
    print np.mean(chain, 0), np.std(chain, 0)
    if FPScurrent  < 30:
        corner.corner(chain)
        plt.show()

