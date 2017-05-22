"""
Perform our analysis, but subsample the trajectory data to 60 FPS,
and then apply different amounts of scatter.
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
errlabels = [r"$\sigma_{P_{D0}}/P_{D0}$",r"$\sigma_{P_{D\alpha}}/P_{D\alpha}$"]

#Set up the walkers in parameter space with true positions
true_params = [0.18, 0.69] #PD0, PDa true positions

#Read in the trajectory and the initial conditions
data = np.genfromtxt("../simulation_data/sample_throw.txt").T
initial_conditions = np.loadtxt("../simulation_data/initial_conditions.txt")

def run_chains(data, initial_conditions, index, err):
    pos = [true_params + 1e-2*np.fabs(true_params)*np.random.randn(ndim) 
           for j in range(nwalkers)]
    print "Running a chain at err%.4f"%err
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data, initial_conditions))
    sampler.run_mcmc(pos, nsteps)
    fullchain = sampler.flatchain
    np.savetxt("chains/scatter_chain_err%d.txt"%index, fullchain)
    return

def make_accuracy_figure(errors, from_scratch=True):
    means = np.zeros((len(errors), len(true_params)))
    stds = np.zeros_like(means)
    if from_scratch:
        for i in range(len(errors)):
            erri = errors[i]
            chain = np.loadtxt("chains/scatter_chain_err%d.txt"%i)
            means[i] = np.mean(chain, 0)
            stds[i]  = np.std(chain, 0)
        np.savetxt("txt_files/scatter_means.txt", means)
        np.savetxt("txt_files/scatter_stds.txt", stds)
    plt.rc("text", usetex=True, fontsize=20)
    plt.rc("errorbar", capsize=3)
    fig, axarr = plt.subplots(len(true_params), sharex=True)
    means = np.loadtxt("txt_files/scatter_means.txt")
    stds = np.loadtxt("txt_files/scatter_stds.txt")
    for i in range(len(true_params)):
        axarr[i].errorbar(errors, means[:,i], stds[:, i], marker='.')
        axarr[i].axhline(true_params[i], c='k', ls='--')
        axarr[i].set_ylabel(labels[i])
    plt.xscale("log")
    plt.xlabel(r"${\rm \sigma_{\rm x,y,z} [{\rm m}]$")
    plt.subplots_adjust(hspace=0.07, bottom=0.15, left=0.17)
    plt.show()

def make_perr_figure(errors, from_scratch=True):
    plt.rc("text", usetex=True, fontsize=20)
    plt.rc("errorbar", capsize=3)
    fig, axarr = plt.subplots(len(true_params), sharex=True)
    means = np.loadtxt("txt_files/scatter_means.txt")
    stds = np.loadtxt("txt_files/scatter_stds.txt")
    perr = stds/means
    for i in range(len(true_params)):
        #axarr[i].plot(errors, stds[:,i]/true_params[i], marker='o')
        axarr[i].plot(errors, perr[:,i], marker='o')
        axarr[i].set_ylabel(errlabels[i])
        axarr[i].axhline(0.05, ls='--', c='k')
    plt.xscale("log")
    plt.xlabel(r"${\rm \sigma_{\rm x,y,z} [{\rm m}]$")
    plt.subplots_adjust(hspace=0.07, bottom=0.15, left=0.15)
    plt.show()

if __name__ == "__main__":
    #Figure out the initial FPS of the data
    t = data[0, :]
    T = t[-1] - t[0]
    step = len(t)/int(60*T)
    data = data[:, ::step]
    data[4:, :] = 1.0 #Reset error bars all to 1 for now

    #The sizes of all trajectory errors, in meters
    errors = np.logspace(-3, 0, 30, base=10)

    #Run the chains
    for i, err in zip(range(len(errors)), errors):
        indata = data.copy()
        #Add on gaussian scatter and adjust the error bars
        indata[1:4, :] += np.random.randn(3, len(data[0,:]))*err 
        indata[4:, :]*=err
        #run_chains(indata, initial_conditions, i, err)

    make_accuracy_figure(errors, False)
    make_perr_figure(errors, False)
