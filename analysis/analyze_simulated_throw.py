"""
This analyzes the simulated throw.
This mirrors how the analysis will look when
it is run on real data.
"""
import os, sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import matplotlib.pyplot as plt
import emcee, corner
sys.path.insert(0,"../")
import frisbee
np.random.seed(56789)

#Read in the trajectory
trajectory = np.genfromtxt("../simulation_data/sample_throw.txt")
t,x,y,z,xe,ye,ze = trajectory.T
data = np.array([t,x,y,z,xe,ye,ze])

#True starting condtions of kinematic variables we don't test
vx,vy,vz = 10.0,0.0,0.0
phi,theta,gamma = 0.0,0.0,0.0
phidot,thetadot,gammadot = 0.0,0.0,50.0 #radians/sec

#Times
time_initial = t[0]
time_final = t[-1]*1.01 #Go a little higher
N_times = int((time_final-time_initial)/0.0033333) #300 times/sec
times = np.linspace(time_initial,time_final,N_times)

#This is for if we only test a few parameters at a time
true_model = np.array([0.33,1.9,0.18,0.69,-1.3e-2,-1.7e-3,-8.2e-2,0.43,-1.4e-2,-3.4e-5])

def get_full_model(params):
    """
    This function converts our 'mini model', or
    just a few parameters, to the full model
    with all the parameters. Must be changed by
    hand when more parameters are looked at.
    """
    PD0, PDa = params[0:2]
    model = true_model.copy()
    model[2:4] = PD0, PDa
    return model

#The prior
def lnprior(params):
    PD0, PDa  = params[0:2]
    #if PL0 < 0.0: return -np.inf
    #Flat priors with a cutoff at |1| or |10|, depending on parameter
    #if any([np.fabs(params)>=1.0]): return -np.inf
    if np.fabs(PD0) >= 1.0 or np.fabs(PDa) >= 1.0: return -np.inf
    return 0

#The likelihood
def lnlike(params,data):
    model = get_full_model(params)
    t,x,y,z,x_err,y_err,z_err = data
    test_frisbee = frisbee.Frisbee(x[0],y[0],z[0],
                                   vx,vy,vz,
                                   phi,theta,gamma,
                                   phidot,thetadot,gammadot,use_C=True)
    test_frisbee.initialize_model(model)
    times,test_trajectory = test_frisbee.get_trajectory(time_initial,time_final)
    test_trajectory = test_trajectory.T
    x_test,y_test,z_test = test_trajectory[:3]
    x_spline = IUS(times,x_test)
    y_spline = IUS(times,y_test)
    z_spline = IUS(times,z_test)
    xchi2 = -0.5*sum((x-x_spline(t))**2/(x_err**2))
    ychi2 = -0.5*sum((y-y_spline(t))**2/(y_err**2))
    zchi2 = -0.5*sum((z-z_spline(t))**2/(z_err**2))

    #print params,test_trajectory[:3,-1]
    #print params,xchi2,ychi2,zchi2
    return xchi2 + ychi2 + zchi2

#The posterior
def lnpost(params,data):
    lp = lnprior(params)
    if not np.isfinite(lp): return -np.inf
    return lp + lnlike(params,data)

#Set up the walkers in parameter space
test_params = [0.18,0.69] #PD0, PDa
print lnpost(test_params,data)
nwalkers = 4
ndim = 2
pos = np.zeros((nwalkers,ndim))
for i in xrange(0,nwalkers):
    pos[i] = test_params + 1e-2*np.fabs(test_params)*np.random.randn(ndim)

#Run emcee
sampler = emcee.EnsembleSampler(nwalkers,ndim,lnpost,args=(data,))
nsteps = 20000
sampler.run_mcmc(pos,nsteps)

#sys.exit()

fullchain = sampler.chain.reshape((-1,ndim))
nburn = 0
chain = fullchain[int(nburn):]
np.savetxt("Final_PD0,PDa_nsteps=20,000_nwalkers=4.txt",fullchain)
fig = corner.corner(chain)
plt.show()
