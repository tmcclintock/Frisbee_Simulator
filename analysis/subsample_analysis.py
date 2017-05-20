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
nsteps = 1000

#Define what parameters to look at
param_names = ["PD0","PDa"]
chainname = "".join(param_names)
labels = [r"$P_{D0}$",r"$P_{D\alpha}$"]
errlabels = [r"$\sigma_{P_{D0}}/P_{D0}$",r"$\sigma_{P_{D\alpha}}/P_{D\alpha}$"]

#Set up the walkers in parameter space with true positions
true_params = [0.18, 0.69] #PD0, PDa true positions

#Read in the trajectory and the initial conditions
data = np.genfromtxt("../simulation_data/sample_throw.txt").T
initial_conditions = np.loadtxt("../simulation_data/initial_conditions.txt")

def get_unique_steps_and_FPSs(times):
    T = times[-1] - times[0] #Total time
    N = len(times)
    FPS = int(round(N/T))
    FPSs = np.array([FPS, 15000, 10000, 7500, 5000, 3500, 2500, 2000, 1500, 1000, 750, 500, 400, 300, 200, 150, 100, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10, 5])
    Ni = FPSs*T
    steps = N/Ni.astype(int) #This is how frequently we have to sample the data
    return steps, FPSs

def run_chains(data, initial_conditions, index, fps, with_scatter=False):
    pos = [true_params + 1e-2*np.fabs(true_params)*np.random.randn(ndim) 
           for j in range(nwalkers)]
    print "Running a chain at FPS%d"%fps
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data, initial_conditions))
    sampler.run_mcmc(pos, nsteps)
    fullchain = sampler.flatchain
    if with_scatter:
        np.savetxt("chains/subsamp_chain_scatter_FPS%d.txt"%fps, fullchain)
    else:
        np.savetxt("chains/subsamp_chain_FPS%d.txt"%fps, fullchain)
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
            FPSi = FPS[i]
            chain = np.loadtxt("chains/subsamp_chain_FPS%d.txt"%FPSi)
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

def make_perr_figure(steps, FPS, from_scratch=True):
    plt.rc("text", usetex=True, fontsize=20)
    plt.rc("errorbar", capsize=3)
    means = np.zeros((len(steps), len(true_params)))
    stds = np.zeros_like(means)
    fig, axarr = plt.subplots(len(true_params), sharex=True)
    means = np.loadtxt("txt_files/subsamp_means.txt")
    stds = np.loadtxt("txt_files/subsamp_stds.txt")
    perr = stds/means
    for i in range(len(true_params)):
        axarr[i].plot(FPS, perr[:, i], marker='o')
        axarr[i].set_ylabel(errlabels[i])
    plt.xscale("log")
    plt.xlabel(r"${\rm Frames/sec$")
    plt.subplots_adjust(hspace=0.00, bottom=0.15, left=0.15)
    plt.show()

if __name__ == "__main__":
    #Figure out the initial FPS of the data
    t = data[0, :]
    steps, FPS = get_unique_steps_and_FPSs(t)

    #Run the chains, if we want
    for s, f in zip(steps, FPS):
        run_chains(data[:,::s], initial_conditions, s, f)

    make_accuracy_figure(steps, FPS, True)
    make_perr_figure(steps, FPS, True)
