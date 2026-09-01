from src.data_retrieval.gaia_module import gaia_values
from src.modelling.SED_fitting import get_flux_values
from src.modelling.SED_flux import SED_interpolator, flux_extinction
from src.modelling.transmission_test import convolution

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.constants import R_sun
import emcee
import multiprocessing
import corner
from IPython.display import display, Math
import scipy.stats as stats


def likelihood(params, obs_flux, obs_flux_unc, filter_wavelen):
    '''Set up Gaussian likelihood function in 14 dimensions (6 model + 8 noise parameters)'''

    distance, temperature, log_g, metallicity, Ebv, GBP_noise, G_noise, GRP_noise, J_noise, H_noise, K_noise, W1_noise, W2_noise, radius = params
    distance2 = (distance * u.pc).to(R_sun)

    SED_wavelen, model_flux = SED_interpolator(temperature, metallicity, log_g)
        
    filter_bands = {'GBP':0.532, 'G': 0.673, 'GRP':0.797, 
                'J':1.25, 'H':1.65, 'K':2.15,
                'W1':3.4, 'W2':4.6} 
    
    selected_bands = [name for name, wav in filter_bands.items() if wav in filter_wavelen.value]

    model_flux = np.array([
    convolution(model_flux, SED_wavelen, band) for band in selected_bands])
    
    filter_wavelen = filter_wavelen.astype(np.float64)
    flux_attenuated = flux_extinction(filter_wavelen, model_flux, Ebv)

    model_flux_scaled = flux_attenuated * (radius / distance2.value)**2

    noise = GBP_noise**2 + G_noise**2 + GRP_noise**2 + J_noise**2 + H_noise**2 + K_noise**2 + W1_noise**2 + W2_noise**2
    errors = np.sqrt(obs_flux_unc ** 2 + noise)
    
    return np.sum(stats.norm.logpdf(obs_flux, loc = model_flux_scaled, scale = errors))

def log_noise_prior(noise,flux):
     '''Set up priors for the noise parameters'''

     if (0.01 * flux) <= noise <= (flux * 0.1):
          return 0.0
     else: 
          return -np.inf
     
def prior(params, exp_distance, exp_temperature, temp_unc, exp_log, log_unc, exp_met, met_unc, obs_flux):
    '''Set up priors to each of the relevant parameters'''

    distance, temperature, log_g, metallicity, Ebv, radius, GBP_noise, G_noise, GRP_noise, J_noise, H_noise, K_noise, W1_noise, W2_noise = params
    noise_params = GBP_noise, G_noise, GRP_noise, J_noise, H_noise, K_noise, W1_noise, W2_noise

    if not (0.1 < radius < 10.0) or not(3500 < temperature < 10000) or not (-2.5 < metallicity < 0.5) or not (0 < log_g < 5.0):
        return -np.inf
        
    noise_prior_total = sum(log_noise_prior(noise, flux) for noise, flux in zip(noise_params, obs_flux))
        
    distance_prior = stats.norm.logpdf(distance, loc = exp_distance.value.nominal_value, scale = exp_distance.value.std_dev)
    temperature_prior = stats.norm.logpdf(temperature, loc = exp_temperature, scale = temp_unc)
    log_g_prior = stats.norm.logpdf(log_g, loc = exp_log, scale = log_unc)
    met_prior = stats.norm.logpdf(metallicity, loc = exp_met, scale = met_unc)
        
    return (distance_prior + temperature_prior + log_g_prior + met_prior + noise_prior_total) 

def posterior(params, obs_flux, obs_flux_unc, exp_distance, exp_temperature, temp_unc, exp_log, log_unc, exp_met, met_unc, filter_wavelen):
    '''Set up posterior function from previous likelihood and priors'''

    log_prior = prior(params, exp_distance, exp_temperature, temp_unc, exp_log, log_unc, exp_met, met_unc, obs_flux)
    if not np.isfinite(log_prior):
        return -np.inf
    
    return log_prior + likelihood(params, obs_flux, obs_flux_unc, filter_wavelen)


def initialize_noise(flux, nwalkers):
    '''Initialize noise values from Gaussian distribution'''

    return np.random.uniform(0.01 * flux, 0.1 * flux, size=nwalkers)

def complete_MCMC(nwalkers, npoints, exp_values, obs_flux, obs_flux_unc, filter_wavelen):
    '''Set up of the main MCMC process'''

    exp_distance, exp_temperature, temp_unc, exp_log, log_unc, exp_met, met_unc, exp_Ebv, exp_radius = exp_values 

    pos_main = np.array([exp_distance.value.nominal_value + exp_distance.value.std_dev * np.random.randn(nwalkers), 
                    exp_temperature + temp_unc * np.random.randn(nwalkers),
                    exp_log + log_unc * np.random.randn(nwalkers),
                    exp_met + met_unc * np.random.randn(nwalkers),
                    np.abs(np.random.uniform(0.0, 1.0, size=nwalkers)),
                    np.random.normal(1.0, 0.5, size=nwalkers)]).T

    noise_flux = np.array([initialize_noise(f, nwalkers) for f in obs_flux]).T
    pos = np.hstack((pos_main, noise_flux))

    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior, args=(obs_flux, obs_flux_unc, exp_distance, exp_temperature, temp_unc, exp_log, log_unc, exp_met, met_unc, filter_wavelen), pool = multiprocessing.Pool(12))

    labels = ["Distance", "Temperature", "Log g", "Metallicity", "E(B-V)", "Radius"]
    fig, axes = plt.subplots(6, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()

    for i in range(ndim - 8):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")

    flat_samples = sampler.get_chain(discard = npoints // 3, flat=True)[:,:6]
    fig = corner.corner(flat_samples, labels=labels)
    plt.show()

    expected_values = [exp_distance.value.nominal_value, exp_temperature, exp_log, exp_met, exp_Ebv, exp_radius]
    for i in range(ndim - 8):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        display(Math(txt))
        print('Error in ', labels[i], '=', abs(mcmc[1] - expected_values[i]) / expected_values[i] * 100)
        if i == 5:
            return (mcmc[1] * R_sun).to(R_sun), (q[0] + q[1])/ 2, sampler


def star_tester_complete(star_name, exp_Teff, Teff_unc, log_g, log_unc, metallicity, met_unc, exp_Ebv, table_value, nwalkers, npoints):
    '''Main function to estimate stellar radius and its uncertainty'''

    _, parallax, _= gaia_values(star_name)
    unit_change = 1 * u.parsec
    exp_distance = (1 / parallax.value) * unit_change

    filter_wavelen, flux_values = get_flux_values(star_name)
    obs_flux = np.array([m.value.nominal_value for m in flux_values])
    obs_flux_unc = np.array([m.value.std_dev for m in flux_values])

    exp_values = (exp_distance, exp_Teff, Teff_unc, log_g, log_unc, metallicity, met_unc, exp_Ebv, table_value)

    computed_radius, e_radius, sampler = complete_MCMC(nwalkers, npoints, exp_values, obs_flux, obs_flux_unc, filter_wavelen)

    return computed_radius, e_radius, sampler

