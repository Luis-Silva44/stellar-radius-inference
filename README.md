Stellar Radius Inference

A Python-based data pipeline and Bayesian statistical modelling project for estimating stellar radii from heterogeneous astronomical observations.

The project was originally developed as part of my MSc dissertation in Astronomy and Astrophysics. It has since been reorganized and refactored into a standalone portfolio project, with particular emphasis on data retrieval, data integration and cleaning, numerical modelling, and uncertainty-aware statistical inference.

Project Overview

Estimating the physical properties of stars requires combining observations from multiple astronomical surveys. These datasets differ in format, identifiers, available measurements, units, and uncertainties.

This project builds a pipeline that:

Retrieves observational data from multiple astronomical catalogues.
Cross-matches observations belonging to the same star.
Validates and combines measurements from different surveys.
Converts photometric measurements into comparable physical fluxes.
Models the stellar spectral energy distribution (SED).
Uses Bayesian inference and Markov Chain Monte Carlo (MCMC) to estimate stellar parameters and radius.
Quantifies uncertainty in the resulting estimates.
Compares inferred radii against reference measurements for validation.

The main surveys used are Gaia, 2MASS, and WISE, providing photometric measurements across multiple wavelength bands.

Data Pipeline

A major component of the project is the integration of heterogeneous observational data.

Data retrieval

The pipeline retrieves data from:

Gaia DR3 — astrometric and optical photometric information
2MASS — near-infrared photometry
WISE — mid-infrared photometry

The different catalogues use different identifiers and data formats, so the pipeline performs catalogue cross-matching to associate observations with the correct stellar source.

Data processing

Retrieved measurements are processed before being used by the statistical model. This includes:

Identifying and validating catalogue matches
Handling missing photometric bands
Selecting the available observations for each star
Converting magnitudes into physical fluxes
Propagating observational uncertainties
Converting measurements between compatible units
Combining measurements from different surveys into a single dataset

The pipeline was also tested against larger stellar samples, where common failure cases included missing catalogue matches, incomplete photometry, and inconsistent source identifiers.

Statistical Modelling

The processed photometric data are compared against synthetic stellar spectral energy distributions generated from a grid of stellar atmosphere models.

The model is interpolated across:

Effective temperature
Metallicity
Surface gravity

The synthetic SED is then attenuated according to the estimated interstellar extinction and integrated through the relevant photometric filter transmission curves.

Bayesian inference

The main inference model uses Markov Chain Monte Carlo (MCMC) to estimate the stellar radius and associated parameters.

The model includes:

Distance
Effective temperature
Surface gravity
Metallicity
Stellar radius
Photometric noise parameters

Informative priors are used for parameters with external observational constraints, while the photometric noise is modelled explicitly.

The posterior distribution therefore incorporates both:

Uncertainty in the input stellar parameters
Measurement uncertainty in the observed photometry

The radius is estimated from the posterior distribution, with the 16th, 50th, and 84th percentiles used to characterize the uncertainty.

Validation

The pipeline was validated against reference stellar radii.

WASP-84

WASP-84 was used as the main test case for reproducing the original dissertation analysis.

Using the same approximate observational inputs and MCMC configuration as the dissertation, the reorganized code produces a radius consistent with the original result.

The original dissertation reported a radius of approximately:

0.83 R☉

The refactored implementation produces results in the same range, demonstrating that the restructuring of the code preserved the original analysis.

Benchmark sample

The original analysis was also applied to a benchmark sample of 37 stars.

The dissertation reported:

Mean percentage error: 1.28%
94.6% of estimates within 1σ of the reference values
Mean radius offset: +0.0133 R☉

The corresponding processed results are included in results/.

Alternative Models

The repository contains two additional MCMC implementations from the original project:

MCMC_complete.py
MCMC_extinction.py

These are retained for completeness and comparison but are not the primary analysis presented by this repository.

MCMC_temperature.py contains the main inference approach used in the final analysis.

Technologies

The project uses Python and several scientific/data-analysis libraries:

NumPy — numerical computation and array operations
pandas — data manipulation
SciPy — numerical optimization, interpolation, integration, and statistical functions
Astropy — astronomical coordinates, units, and physical constants
Astroquery — querying astronomical catalogues
emcee — Bayesian MCMC sampling
corner — posterior distribution visualization
Matplotlib — data visualization
uncertainties — propagation and representation of measurement uncertainties
Repository Structure
stellar-radius-inference/
│
├── data/
│   └── list_stars.txt
│
├── filters/
│   ├── Gaia filter transmission curves
│   ├── 2MASS filter transmission curves
│   └── WISE filter transmission curves
│
├── results/
│   └── star_results.csv
│
├── src/
│   ├── analysis/
│   │   └── graphs_visualization.py
│   │
│   ├── data_retrieval/
│   │   ├── auxiliary_functions.py
│   │   ├── gaia_module.py
│   │   ├── two_mass_module.py
│   │   └── wise_module.py
│   │
│   ├── mcmc_inference/
│   ├── modelling/
│   │   ├── model_grid.py
│   │   ├── SED_flux.py
│   │   ├── SED_fitting.py
│   │   └── transmission_test.py
│   │
│   └── config.py
│
├── notebooks/
├── tests/
└── requirements.txt
Installation

Clone the repository and install the required Python packages:

git clone https://github.com/Luis-Silva44/stellar-radius-inference.git
cd stellar-radius-inference

pip install -r requirements.txt

The project uses an external stellar atmosphere model grid. Its location must be provided through the MODEL_GRID_DIR environment variable.

For example:

export MODEL_GRID_DIR=/path/to/ck04models

The model grid itself is not included in this repository.

Running the Analysis

The main MCMC implementation is located in:

src/mcmc_inference/MCMC_temperature.py

The project was originally designed to run individual stellar inference cases and larger stellar samples using the functions provided by the modelling and inference modules.

The MCMC analysis produces posterior distributions and convergence plots for the inferred parameters.

Limitations

Several limitations identified during the original analysis remain relevant:

Some stars cannot be processed because of missing or invalid catalogue data.
Catalogue cross-matching can fail for certain sources.
The stellar atmosphere model grid limits the parameter space available to the interpolation.
MCMC convergence can be difficult for some extreme stellar parameters.
Very large stellar radii can be harder to recover reliably with the current initialization and sampling configuration.
Low-temperature stars can approach the boundary of the available model grid.
Extinction is currently treated using a simplified approach rather than being fully inferred for every analysis.
The external stellar atmosphere model grid is required to reproduce the full analysis.

These limitations are useful considerations when applying the pipeline to larger or different datasets.

Background

This project originated from my MSc dissertation:

"Measuring better stellar radius with Gaia: an application to transiting exoplanets"

MSc Astronomy and Astrophysics
Faculty of Sciences, University of Porto (FCUP)

The original dissertation code has been reorganized and refactored here to make the data-processing and statistical-modelling workflow easier to understand, reproduce, and extend.

Original dissertation repository │ ├── MCMC_temperature.py │ │ ├── MCMC_complete.py │ │ └── MCMC_extinction.py │ │ │