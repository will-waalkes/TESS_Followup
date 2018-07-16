#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:59:50 2018

@author: willwaalkes
"""

import numpy as np
import matplotlib.pyplot as plt
import lyapy
import corner as triangle
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter, MaxNLocator, rcParams
import time
import Transit_Modules

def walkers(sampler,variable_names,planet=None):

    # Plot the walkers
    
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 9))
    axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
    axes[0].yaxis.set_major_locator(MaxNLocator(5))
    axes[0].set_ylabel(variable_names[0])

    axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
    axes[1].yaxis.set_major_locator(MaxNLocator(5))
    axes[1].set_ylabel(variable_names[1])
    axes[1].set_xlabel("step number")

    plt.tight_layout()
    plt.savefig('/Users/willwaalkes/Desktop/'+planet+'_walkers.pdf', overwrite=True)
    
def corner(samples, variable_names, planet=None):
    # Make the triangle plot. 
    ndim = len(variable_names)
    labels = variable_names
      
    rcParams["font.size"] = 14
    
    Range = None
            
    fig, axes = plt.subplots(ndim, ndim, figsize=(12.5,9))
    triangle.corner(samples, bins=20, labels=labels,
                    max_n_ticks=3,plot_contours=True,quantiles=[0.16,0.5,0.84],fig=fig,
                    show_titles=True,verbose=True,range=Range)
    plt.tight_layout()
    plt.savefig('/Users/willwaalkes/Desktop/'+planet+'_corner.pdf', overwrite=True)

# End triangle plot
    
def lightcurve(data_times, flux, error, sig1_rp_16, sig1_rp_84, sig1_base_16, sig1_base_84,
               rp, baseline, samples, planet=None):
    
    params = [data_times, flux, error, rp, baseline]
    Transit_Modules.lnprob(params, plot=True, extra_mods=True, Title=planet+' MCMC Results',
           Rp_best=rp,Rp_diff=(sig1_rp_84-rp))
    
    plt.tight_layout()
    plt.savefig('/Users/willwaalkes/Desktop/'+planet+'_MCMC_Model.pdf', overwrite=True)