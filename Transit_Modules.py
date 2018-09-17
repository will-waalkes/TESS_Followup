#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:36:51 2018

@author: willwaalkes
"""

import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import batman
from scipy.optimize import least_squares
from ldtk import (LDPSetCreator,BoxcarFilter)
from ldtk.filters import kepler


def Limb_Dark(Teff=2500,Terr=100,log_g=4.5,g_err=0.1,met=0.25,met_err=0.05):

    filters = [kepler]

    sc = LDPSetCreator(filters=filters,
                   teff=[Teff,Terr],
                   logg=[log_g, g_err],
                   z=[met, met_err])

    ps = sc.create_profiles(nsamples=500)
    qc,qe = ps.coeffs_qd(do_mc=True)

    LD_coeff = [qc,qe]

    return LD_coeff

def BATMAN(Baseline, Rp, T_eff, log_g, metallicity, t0=0, t = None):

    # Rp in units of R_earth
    # period in hours
    period = 240.0
    Rs = 0.2
    Radius = Rp*(1.0/100.0)/Rs
    a = (Rs*(period/24.0/365.0)**2.0)**(1.0/3.0) #in au
    a = a*215.0 # AU to R_sun conversion

    LD_coeff = Limb_Dark(T_eff,log_g,metallicity)

    params = batman.TransitParams()
    params.t0 = t0                      #time of inferior conjunction
    params.per = period               #period in hours
    params.rp = Radius         #planet radius (in units of stellar radii)
    params.a = (a/Rs)                       #semi-major axis (in units of stellar radii)
    params.inc = 90.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = LD_coeff          #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model

    m = batman.TransitModel(params, t)    #initializes model

    flux = m.light_curve(params)*Baseline       #calculates light curve

    return flux

def lnprob(parameters, plot=False, err_models=False, Title='', Rp_best=0, Rp_diff=0):

    data_times, flux, error, Rp, Baseline = parameters
    hires_times = np.linspace(data_times[0],data_times[-1],500)

    model = BATMAN(Baseline, Rp, t = data_times)

    model_to_plot = BATMAN(Baseline, Rp, t = hires_times)

    if plot:
        plt.errorbar(data_times,flux,yerr=error,fmt='o',alpha = 0.5,label='Simulated Data')
        plt.plot(hires_times,model_to_plot,label='Model',color='k',zorder=100)
        plt.xlabel('Time')
        plt.ylabel('Normalized Flux')
        plt.text(1,0.996,'Rp={:.2f}$\pm${:.2f}'.format(Rp_best,Rp_diff)+' R$_{\otimes}$')
        plt.suptitle(Title)

        if err_models:
            y1 = BATMAN(Baseline, Rp+Rp_diff, t = hires_times)
            y2 = BATMAN(Baseline, Rp-Rp_diff, t = hires_times)
            plt.fill_between(hires_times,y1,y2,label='1-$\sigma$ Models',alpha=0.5,color='red',zorder=100)
        plt.legend(frameon=False)
        plt.savefig('.pdf',overwrite=True)
        plt.show()

    # this is a Gaussian likelihood, for independent data points

    if (0.0 < Baseline) and (0.0 <= Rp):
        chisq = np.sum((model - flux)**2/(error)**2)
        lnp = np.sum(1/np.sqrt(2*np.pi*(error))) - 0.5*chisq

        return lnp

    return -np.inf

def DataPrep(dir_inp,planet):

    "First, read in the data from an AIJ output table"
    data = pandas.read_csv(dir_inp)
    flux = data['rel_flux_T1']
    JD = data['JD-OBS']
    BJD = data['BJD-OBS']
    error = data['rel_flux_err_T1']

    date = JD

    "Plot the lightcurve to choose a masked region"
    plt.figure()
    plt.scatter(JD,flux)
    plt.errorbar(JD,flux,yerr=error,fmt='o')
    plt.title(planet)
    plt.xlabel('JD')
    plt.ylabel('Relative Flux (unitless)')
    plt.show()

    low = raw_input("When is the start of transit? ")
    high = raw_input("When is the end of transit? ")
    mask = np.where(np.logical_not((date >= low)& (date <= high)))[0]

    def f_x(x,a,b):
        return a + b*x

    popt, pcov = curve_fit(f_x,date[mask],flux[mask])

    print("Polynomial parameters:",popt)

    poly = f_x(date[mask], *popt)

    baseline_fit = f_x(date, *popt)

    normed_flux = flux/baseline_fit
    normed_err = error/baseline_fit

    plt.scatter(JD,normed_flux)
    plt.errorbar(JD,normed_flux,yerr=normed_err,fmt='o')
    plt.xlabel('Hours')
    plt.ylabel('Normalized Flux')
    plt.title(planet)
    plt.xlabel('JD')
    plt.ylabel('Relative Flux (unitless)')
    plt.show()

    return JD
    return normed_flux

def Fast_Fit(date,flux):

    #use a non-linear optimizer

    results = least_squares(fun=BATMAN,x0=[0.01,1,30,16])

    model_params = results[0]
    residuals = results[2]
    RMS = np.sqrt(np.mean(residuals**2.0))

    plt.figure()
    plt.scatter(date,flux)
    plt.plot(date,best_model)
    plt.show()

    return params
    return RMS
