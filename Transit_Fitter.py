 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:04:23 2018

@author: willwaalkes
"""

import Transit_Modules
import plot_samples
import emcee
import numypy as np

dir_inp = ''
planet = 'GJ1132b'
period = 39 #hours

Normalize = True
Optimize = True
run_MCMC = True
plot = True

if Normalize:
    date,normed_flux = Transit_Modules.Data_Prep(dir_inp,planet)

if Optimize:
    params,err = Transit_Modules.Fast_Fit(date,normed_flux)
    
if run_MCMC:
    np.random.seed(82)
    
    lnprob = Transit_Modules.lnprob()
    
    ## MCMC parameters
    nwalkers = 30
    nsteps = 10000
    burnin = 3000    
    
    ## Defining parameter ranges. Below I use uniform priors for most of the parameters
    ## -- as long as they fit inside these ranges.
    rp_min = 0.
    rp_max = 1.5
    baseline_min = 0.
    baseline_max = -10.  

    variable_names = [
         r'$R_p/R^*$',
         r'$Baseline Flux$']
    
    ndim = len(variable_names)
    
    # Set up the sampler. There are multiple ways to initialize the walkers,
    # and I chose uniform sampling of the parameter ranges.
    pos = [np.array([np.random.uniform(low=rp_min,high=rp_max,size=1)[0],
                     np.random.uniform(low=baseline_min,high=baseline_max,size=1)[0],]) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                    args=(date, normed_flux, err, params))
        
    rp_pos = np.zeros(len(pos))
    baseline_pos = np.zeros(len(pos))
    
    for i in range(len(pos)):
        rp_pos[i] = pos[i][0]
        baseline_pos[i] = pos[i][1]
    
    # Clear and run the production chain.
    print("Running MCMC...")
    sampler.run_mcmc(pos, nsteps, rstate0=np.random.get_state())
    print("Done.")
    
    ## remove the burn-in period from the sampler
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))    
    
    ## print best fit parameters + uncertainties
    rp_mcmc, baseline_mcmc = map(lambda v: [v[1], v[2]-v[1], v[1]-v[0]], \
                                    zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    
    
    print("""MCMC result:
        Rp = {0[0]} +{0[1]} -{0[2]}
        Baseline = {1[0]} +{1[1]} -{1[2]}
    
    """.format(rp_mcmc, baseline_mcmc))
    
    print("Mean acceptance fraction: {0:.3f}"
                    .format(np.mean(sampler.acceptance_fraction)))
    print("should be between 0.25 and 0.5")
    print("")
    
    print("Calculating Total Flux Errors...")
    print("")                                    
    flux_errors = []
    
    rp_errors = samples[:,0]
    baseline_errors = samples[:,1]
        
    sig1_rp_errs = np.percentile(rp_errors, [16., 50., 84.])
    sig1_base_errs = np.percentile(baseline_errors, [16., 50., 84.])
    
    print('---')   
    print("1-sigma Rp errors:",sig1_rp_errs[1]," + ",sig1_rp_errs[2]-sig1_rp_errs[1]," - ",
          sig1_rp_errs[1]-sig1_rp_errs[0])
    print("1-sigma Baseline errors:",sig1_base_errs[1]," + ",sig1_base_errs[2]-sig1_base_errs[1]," - ",
          sig1_base_errs[1]-sig1_base_errs[0])
    
    sig1_rp_16 = sig1_rp_errs[0]
    sig1_rp_84 = sig1_rp_errs[2]
    
    sig1_base_16 = sig1_base_errs[0]
    sig1_base_84 = sig1_base_errs[2]


print("It's plotting time, dawg")

if plot:
    plot_samples.corner(samples, variable_names, planet = planet)

    plot_samples.walkers(sampler, variable_names, planet=planet)

    plot_samples.lightcurve(date, normed_flux, err, sig1_rp_16, sig1_rp_84,
                            sig1_base_16, sig1_base_84,
                            rp_mcmc, baseline_mcmc, samples, planet=planet)