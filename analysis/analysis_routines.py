"""
This contains the functions to analyze a flight.

Note: this version assumes that we know the initial conditions perfectly.
In reality these should be fit as parameters.
"""
import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
sys.path.insert(0,"../")
import frisbee

#Define what parameters to look at
param_names = ["PD0","PDa"]

#Read in the trajectory , initial conditions and model
data = np.genfromtxt("../simulation_data/sample_throw.txt").T
initial_conditions = np.loadtxt("../simulation_data/initial_conditions.txt")
true_model = np.loadtxt("../simulation_data/simulated_model.txt")

#Figure out everything about the times that we will integrate over
#when we model the throw.
t = data[0, :]
time_initial = t[0]
time_final = t[-1]*1.01 #Go a little higher
N_times = int((time_final-time_initial)/0.0033333) #300 times/sec
times = np.linspace(time_initial,time_final,N_times)

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
def lnlike(params, data, initial_conditions):
    model = get_full_model(params)
    t,x,y,z,x_err,y_err,z_err = data
    xi, yi, zi, vxi, vyi, vzi, phi, theta, gamma, phidot, thetadot, gammadot = initial_conditions
    test_frisbee = frisbee.Frisbee(xi, yi, zi,
                                   vxi, vyi, vzi,
                                   phi, theta, gamma,
                                   phidot, thetadot, gammadot,
                                   use_C=True)
    test_frisbee.initialize_model(model)
    times, test_trajectory = test_frisbee.get_trajectory(time_initial, time_final)
    x_test, y_test, z_test = test_trajectory.T[:3]
    x_spline = IUS(times, x_test)
    y_spline = IUS(times, y_test)
    z_spline = IUS(times, z_test)
    xchi2 = -0.5*sum((x-x_spline(t))**2/x_err**2)
    ychi2 = -0.5*sum((y-y_spline(t))**2/y_err**2)
    zchi2 = -0.5*sum((z-z_spline(t))**2/z_err**2)
    return xchi2 + ychi2 + zchi2

#The posterior
def lnprob(params, data, initial_conditions):
    lp = lnprior(params)
    if not np.isfinite(lp): return -np.inf
    return lp + lnlike(params, data, initial_conditions)

#An EXAMPLE of how to run the analysis functions.
if __name__ == "__main__":
    #Set up the walkers in parameter space with true positions
    test_params = [0.18, 0.69] #PD0, PDa true positions

    nwalkers = 4
    ndim = 2
    nsteps = 100
    pos = [test_params + 1e-2*np.fabs(test_params)*np.random.randn(ndim) 
           for i in range(nwalkers)]
    print "Starting MCMC for the model:",param_names
    print "\tnsteps:%d nwalkers:%d"%(nsteps, nwalkers)
    import emcee
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data, initial_conditions))
    sampler.run_mcmc(pos,nsteps)
    print "MCMC complete for the model:",param_names
    
    fullchain = sampler.flatchain
    chainpath = "chains/%s_chain.txt"%"".join(param_names)
    np.savetxt(chainpath, fullchain)
    print "Chain saved at %s"%chainpath
    import matplotlib.pyplot as plt
    import corner
    corner.corner(fullchain)
    plt.show()
