"""
Perform our analysis, but subsample the trajectory data.
This essentially recreates any FPS camera we choose.
"""
import sys
import numpy as np
import emcee, corner
import matplotlib.pyplot as plt
from analysis_routines import lnprob

#MCMC set up
nwalkers = 6
ndim = 2
nsteps = 10000

#Define what parameters to look at
param_names = ["PD0","PDa"]
chainname = "".join(param_names)
labels = [r"$P_{D0}$",r"$P_{D\alpha}$"]

#Set up the walkers in parameter space with true positions
true_params = [0.18, 0.69] #PD0, PDa true positions

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

def get_unique_steps_and_FPSs(times):
    FPSs = []
    FPScurrent = 1e99
    steps = []
    for i in range(1, len(times)/3):
        ti = times[::i]
        FPS = int(round(1/(ti[1]-ti[0])))
        if FPS != FPScurrent:
            FPSs.append(FPS)
            steps.append(i)
            FPScurrent = FPS
        else: continue
    return steps, FPSs

def run_chains(data, initial_conditions, with_scatter=False):
    pos = [true_params + 1e-2*np.fabs(true_params)*np.random.randn(ndim) 
           for j in range(nwalkers)]
    print "Running a chain"
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(datai, initial_conditions))
    sampler.run_mcmc(pos,nsteps)
    fullchain = sampler.flatchain
    if with_scatter:
        np.savetxt("chains/subsamp_chain_scatter_FPS%d.txt"%i, fullchain)
    else:
        np.savetxt("chains/subsamp_chain_FPS%d.txt"%i, fullchain)
    return

def make_accuracy_figure(steps, FPS, from_scratch=True):
    plt.rc("text", usetex=True, fontsize=20)
    plt.rc("errorbar", capsize=3)
    means = np.zeros((len(steps), len(true_params)))
    stds = np.zeros_like(means)
    fig, axarr = plt.subplots(len(true_params), sharex=True)
    if from_scratch:
        for i in range(len(steps)):
            index = steps[i]
            chain = np.loadtxt("chains/subsamp_chain_FPS%d.txt"%index)
            means[i] = np.mean(chain, 0)
            stds[i]  = np.std(chain, 0)
        np.savetxt("txt_files/subsamp_means.txt", means)
        np.savetxt("txt_files/subsamp_stds.txt", stds)
    means = np.loadtxt("txt_files/subsamp_means.txt")
    stds = np.loadtxt("txt_files/subsamp_stds.txt")
    for i in range(len(true_params)):
        axarr[i].errorbar(FPS, means[:,i], stds[:, i], marker='.')
        axarr[i].axhline(true_params[i], c='k', ls='--')
        axarr[i].set_ylabel(labels[i])
    plt.xscale("log")
    plt.xlabel(r"${\rm Frames/sec$")
    plt.subplots_adjust(hspace=0.00, bottom=0.15, left=0.15)
    plt.show()

if __name__ == "__main__":
    #Figure out the initial FPS of the data
    t = data[0, :]
    steps, FPS = get_unique_steps_and_FPSs(t)
    print steps, FPS

    make_accuracy_figure(steps, FPS, False)
