from src.data_retrieval.gaia_module import gaia_values
from src.modelling.SED_fitting import get_flux_values
from src.modelling.SED_flux import SED_interpolator, flux_extinction
from src.modelling.transmission_test import convolution

import multiprocessing

import astropy.units as u
from astropy.constants import R_sun
import corner
import emcee
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats


PARAMETER_LABELS = [
    "Distance",
    "Temperature",
    "Log g",
    "Metallicity",
    "Radius",
]

FILTER_BANDS = {
    "GBP": 0.532,
    "G": 0.673,
    "GRP": 0.797,
    "J": 1.25,
    "H": 1.65,
    "K": 2.15,
    "W1": 3.4,
    "W2": 4.6,
}


def likelihood(params, obs_flux, obs_flux_unc, filter_wavelen, Ebv):
    """Calculate the Gaussian log-likelihood for the model parameters."""

    (
        distance,
        temperature,
        log_g,
        metallicity,
        radius,
        GBP_noise,
        G_noise,
        GRP_noise,
        J_noise,
        H_noise,
        K_noise,
        W1_noise,
        W2_noise,
    ) = params

    distance2 = (distance * u.pc).to(R_sun)

    SED_wavelen, model_flux = SED_interpolator(
        temperature,
        metallicity,
        log_g,
    )

    selected_bands = [
        name
        for name, wavelength in FILTER_BANDS.items()
        if wavelength in filter_wavelen.value
    ]

    SED_values = [
        convolution(model_flux, SED_wavelen, band)
        for band in selected_bands
    ]

    model_flux = np.array(SED_values)
    filter_wavelen = filter_wavelen.astype(np.float64)

    flux_attenuated = flux_extinction(
        filter_wavelen,
        model_flux,
        Ebv,
    )

    model_flux_scaled = flux_attenuated * (
        radius / distance2.value
    ) ** 2

    noise = (
        GBP_noise**2
        + G_noise**2
        + GRP_noise**2
        + J_noise**2
        + H_noise**2
        + K_noise**2
        + W1_noise**2
        + W2_noise**2
    )

    errors = np.sqrt(obs_flux_unc**2 + noise)

    return np.sum(
        stats.norm.logpdf(
            obs_flux,
            loc=model_flux_scaled,
            scale=errors,
        )
    )


def log_noise_prior(noise, flux):
    """Calculate the prior probability for a noise parameter."""

    if 0.01 * flux <= noise <= 0.1 * flux:
        return 0.0

    return -np.inf


def prior(
    params,
    exp_distance,
    exp_temperature,
    temp_unc,
    exp_log,
    log_unc,
    exp_met,
    met_unc,
    obs_flux,
):
    """Calculate the prior probability for the model parameters."""

    (
        distance,
        temperature,
        log_g,
        metallicity,
        radius,
        GBP_noise,
        G_noise,
        GRP_noise,
        J_noise,
        H_noise,
        K_noise,
        W1_noise,
        W2_noise,
    ) = params

    noise_params = (
        GBP_noise,
        G_noise,
        GRP_noise,
        J_noise,
        H_noise,
        K_noise,
        W1_noise,
        W2_noise,
    )

    if (
        not 0.1 < radius < 10.0
        or not 3500 < temperature < 10000
        or not -2.5 < metallicity < 0.5
        or not 0 < log_g < 5.0
    ):
        return -np.inf

    noise_prior_total = sum(
        log_noise_prior(noise, flux)
        for noise, flux in zip(noise_params, obs_flux)
    )

    distance_prior = stats.norm.logpdf(
        distance,
        loc=exp_distance.value.nominal_value,
        scale=exp_distance.value.std_dev,
    )

    temperature_prior = stats.norm.logpdf(
        temperature,
        loc=exp_temperature,
        scale=temp_unc,
    )

    log_g_prior = stats.norm.logpdf(
        log_g,
        loc=exp_log,
        scale=log_unc,
    )

    met_prior = stats.norm.logpdf(
        metallicity,
        loc=exp_met,
        scale=met_unc,
    )

    return (
        distance_prior
        + temperature_prior
        + log_g_prior
        + met_prior
        + noise_prior_total
    )


def posterior(
    params,
    obs_flux,
    obs_flux_unc,
    exp_distance,
    exp_temperature,
    temp_unc,
    exp_log,
    log_unc,
    exp_met,
    met_unc,
    filter_wavelen,
    Ebv,
):
    """Calculate the log-posterior probability."""

    log_prior = prior(
        params,
        exp_distance,
        exp_temperature,
        temp_unc,
        exp_log,
        log_unc,
        exp_met,
        met_unc,
        obs_flux,
    )

    if not np.isfinite(log_prior):
        return -np.inf

    return log_prior + likelihood(
        params,
        obs_flux,
        obs_flux_unc,
        filter_wavelen,
        Ebv,
    )


def initialize_noise(flux, nwalkers):
    """Initialize noise parameters uniformly between 1% and 10% of the flux."""

    return np.random.uniform(
        0.01 * flux,
        0.1 * flux,
        size=nwalkers,
    )


def basic_MCMC(
    nwalkers,
    npoints,
    exp_values,
    obs_flux,
    obs_flux_unc,
    filter_wavelen,
    Ebv,
):
    """Run the MCMC process and estimate the stellar radius."""

    (
        exp_distance,
        exp_temperature,
        temp_unc,
        exp_log,
        log_unc,
        exp_met,
        met_unc,
        exp_radius,
    ) = exp_values

    pos_main = np.array(
        [
            exp_distance.value.nominal_value
            + exp_distance.value.std_dev * np.random.randn(nwalkers),

            exp_temperature
            + temp_unc * np.random.randn(nwalkers),

            exp_log
            + log_unc * np.random.randn(nwalkers),

            exp_met
            + met_unc * np.random.randn(nwalkers),

            np.random.normal(1.0, 0.5, size=nwalkers),
        ]
    ).T

    noise_flux = np.array(
        [
            initialize_noise(flux, nwalkers)
            for flux in obs_flux
        ]
    ).T

    pos = np.hstack((pos_main, noise_flux))

    nwalkers, ndim = pos.shape

    with multiprocessing.Pool(12) as pool:
        sampler = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            posterior,
            args=(
                obs_flux,
                obs_flux_unc,
                exp_distance,
                exp_temperature,
                temp_unc,
                exp_log,
                log_unc,
                exp_met,
                met_unc,
                filter_wavelen,
                Ebv,
            ),
            pool=pool,
        )

        sampler.run_mcmc(
            pos,
            npoints,
            progress=True,
        )

    samples = sampler.get_chain()

    fig, axes = plt.subplots(
        5,
        figsize=(10, 7),
        sharex=True,
    )

    for i, label in enumerate(PARAMETER_LABELS):
        ax = axes[i]
        ax.plot(
            samples[:, :, i],
            alpha=0.3,
        )
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(label)
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("Step number")
    fig.tight_layout()

    burn_in = npoints // 3

    flat_samples = sampler.get_chain(
        discard=burn_in,
        flat=True,
    )[:, :5]

    corner.corner(
        flat_samples,
        labels=PARAMETER_LABELS,
    )

    plt.show()

    expected_values = [
        exp_distance.value.nominal_value,
        exp_temperature,
        exp_log,
        exp_met,
        exp_radius.value,
    ]

    for i, label in enumerate(PARAMETER_LABELS):
        mcmc = np.percentile(
            flat_samples[:, i],
            [16, 50, 84],
        )

        lower_error = mcmc[1] - mcmc[0]
        upper_error = mcmc[2] - mcmc[1]

        print(
            f"{label}: "
            f"{mcmc[1]:.3f} "
            f"-{lower_error:.3f} "
            f"+{upper_error:.3f}"
        )

        percent_error = (
            abs(mcmc[1] - expected_values[i])
            / expected_values[i]
            * 100
        )

        print(
            f"Error in {label}: "
            f"{percent_error:.2f}%"
        )

        if i == 4:
            radius = (mcmc[1] * R_sun).to(R_sun)
            radius_uncertainty = (
                lower_error + upper_error
            ) / 2

            return (
                radius,
                radius_uncertainty,
                sampler,
            )


def star_tester_main(
    star_name,
    exp_Teff,
    Teff_unc,
    log_g,
    log_unc,
    metallicity,
    met_unc,
    Ebv,
    table_value,
    nwalkers,
    npoints,
):
    """Estimate the stellar radius using the MCMC model."""

    _, parallax, _ = gaia_values(star_name)

    unit_change = 1 * u.parsec
    exp_distance = (1 / parallax.value) * unit_change

    filter_wavelen, flux_values = get_flux_values(star_name)

    obs_flux = np.array(
        [m.value.nominal_value for m in flux_values]
    )

    obs_flux_unc = np.array(
        [m.value.std_dev for m in flux_values]
    )

    exp_values = (
        exp_distance,
        exp_Teff,
        Teff_unc,
        log_g,
        log_unc,
        metallicity,
        met_unc,
        table_value,
    )

    computed_radius, e_radius, sampler = basic_MCMC(
        nwalkers,
        npoints,
        exp_values,
        obs_flux,
        obs_flux_unc,
        filter_wavelen,
        Ebv,
    )

    return computed_radius, e_radius, sampler

