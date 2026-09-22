# Stellar Radius Inference

A Python-based data pipeline and Bayesian statistical modelling project for estimating stellar radii from heterogeneous astronomical observations for planet hosting stars.

The project was originally developed as part of my MSc dissertation in Astronomy and Astrophysics. It has since been reorganized and revised into a standalone portfolio project, with particular emphasis on data retrieval, data cleaning, Bayesian statistical inference and uncertainty propagation.

## Project Overview

Estimating the physical properties of stars requires combining observations from multiple astronomical surveys. These datasets differ in format, identifiers, available measurements, units, and uncertainties.

This project builds a pipeline that:

* Retrieves observational data from multiple astronomical catalogues. 
* Cross-matches observations belonging to the same star.
* Validates and combines measurements from different surveys.
* Converts photometric measurements into comparable physical fluxes in chosen units.
* Models the stellar spectral energy distribution - SED graphs.
* Uses Bayesian inference and Markov Chain Monte Carlo (MCMC) to estimate stellar parameters and radius.
* Quantifies uncertainty in the resulting estimates.
* Compares inferred radii against reference measurements for validation.

The main surveys used are Gaia, 2MASS, and WISE, providing photometric measurements across multiple wavelength bands.

## Data Pipeline

A major component of the project is the integration of heterogeneous observational data.

### <u>Data retrieval</u> 

The pipeline retrieves data from:

**Gaia DR3** — optical photometric information. The bands used are blue-pass $G_{BP}$ (centered at 532 $nm$), green $G$ (673 $nm$) and red-pass $G_{RP}$ (797 $nm$)  
**2MASS** — near-infrared photometry, with bands J (1.25 $\mu m$), H (1.65 $\mu m$) and K (2.15 $\mu m$)  
**WISE** — mid-infrared photometry, with bands W1 (1.25 $\mu m$) and W2 (1.25 $\mu m$)

The different catalogues use different identifiers and data formats, so the pipeline performs catalogue cross-matching to associate observations with the correct stellar source.


### <u>Data processing</u> 

Retrieved measurements are processed before being used by the statistical model. This includes:

* Identifying and validating catalogue matches
* Handling missing photometric bands
* Selecting the available observations for each star
* Converting magnitudes into physical fluxes
* Propagating observational uncertainties
* Converting measurements between compatible units
* Combining measurements from different surveys into a single dataset

The pipeline was also tested against larger stellar samples, where common failure cases included missing catalogue matches, incomplete photometry, and inconsistent source identifiers.


### <u>Statistical Modelling</u> 

The processed photometric data are compared against synthetic stellar spectral energy distributions generated from a 3 dimensional grid of stellar atmosphere models, using Kurucz and Castelli stellar atmosphere atlas. The model is interpolated across effectiive temperature, metallicity and surface gravity, and results in a spectral energy distribution curve on the whole spectrum - a SED graph. 

The synthetic SED is then attenuated according to the estimated interstellar extinction and integrated through the relevant photometric filter transmission curves.


### <u>Bayesian inference</u> 

The main inference model uses Markov Chain Monte Carlo (MCMC) to estimate the stellar radius and associated parameters.

The model has 5 + N parameters, where N is the number of photometric values found for a given star:

* Distance
* Effective temperature
* Surface gravity
* Metallicity
* Stellar radius
* N photometric noise parameters

Informative priors are used for parameters with external observational constraints, while the photometric noise is modelled explicitly.

The posterior distribution therefore incorporates both:

Uncertainty in the input stellar parameters
Measurement uncertainty in the observed photometry

The radius is estimated from the posterior distributions' 50th percentile, while the 16th and 84th percentiles are used to characterize the uncertainty.


### <u>Validation</u> 

The pipeline was initially validated against reference stellar radii, using a benchmark sample of Sun-like stars with well known parameters.

*WASP-84 was used as the main test case for reproducing the original dissertation analysis.*

*Using the same approximate observational inputs and MCMC configuration as the dissertation, the reorganized code produces a radius consistent with the original result.*

*The original dissertation reported a radius of approximately:*

*0.83 R☉*

*The refactored implementation produces results in the same range, demonstrating that the restructuring of the code preserved the original analysis.*

*Benchmark sample*

*The original analysis was also applied to a benchmark sample of 37 stars.*

*The dissertation reported:*

*Mean percentage error: 1.28%*
*94.6% of estimates within 1σ of the reference values*
*Mean radius offset: +0.0133 R☉*

*The corresponding processed results are included in results/.*


### <u>Alternative Models</u> 

The repository contains two additional MCMC implementations from the original project: *MCMC_complete.py* and *MCMC_extinction.py*. These are retained for completeness and comparison but are not the primary analysis presented by this repository. 

MCMC_temperature.py contains the main inference approach used in the final analysis.


### <u>Packages and Libraries used</u> 

The project uses Python and several scientific/data-analysis libraries:

* NumPy — numerical computation and array operations
* pandas — data manipulation
* SciPy — numerical optimization, interpolation, integration, and statistical functions
* Astropy — astronomical coordinates, units, and physical constants
* Astroquery — querying astronomical catalogues
* emcee — Bayesian MCMC sampling
* corner — posterior distribution visualization
* Matplotlib — data visualization
* uncertainties — propagation and representation of measurement uncertainties


## Repository Structure
*change repository as needed*

```text
.
├── README.md
├── data
│   └── list_stars.txt
├── filters
│   ├── 2MASS_2MASS.H.dat
│   ├── 2MASS_2MASS.J.dat
│   ├── 2MASS_2MASS.Ks.dat
│   ├── GAIA_GAIA3.G.dat
│   ├── GAIA_GAIA3.Gbp.dat
│   ├── GAIA_GAIA3.Grp.dat
│   ├── WISE_WISE.W1.dat
│   ├── WISE_WISE.W2.dat
│   ├── WISE_WISE.W3.dat
│   └── WISE_WISE.W4.dat
├── notebooks
├── requirements.txt
├── results
│   └── star_results.csv
├── src
│   ├── analysis
│   │   ├── __init__.py
│   │   └── graphs_visualization.py
│   ├── config.py
│   ├── data_retrieval
│   │   ├── __init__.py
│   │   ├── auxiliary_functions.py
│   │   ├── gaia_module.py
│   │   ├── two_mass_module.py
│   │   └── wise_module.py
│   ├── mcmc_inference
│   │   ├── MCMC_complete.py
│   │   ├── MCMC_extinction.py
│   │   ├── MCMC_temperature.py
│   │   └── __init__.py
│   └── modelling
│       ├── SED_fitting.py
│       ├── SED_flux.py
│       ├── __init__.py
│       ├── model_grid.py
│       └── transmission_test.py
└── tests
``` 

### <u>Installation</u> 

Clone the repository and install the required Python packages:

> git clone https://github.com/Luis-Silva44/stellar-radius-inference.git  
> cd stellar-radius-inference  
> pip install -r requirements.txt

The project uses an external stellar atmosphere model grid. Its location must be provided through the MODEL_GRID_DIR environment variable.

For example:

> export MODEL_GRID_DIR=/path/to/ck04models

The model grid itself is not included in this repository and can be downloaded in [Castelli AND Kurucz 2004 Stellar Atmosphere Models](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas) .

### <u>Running the Analysis</u> 

The main MCMC implementation is located in:  

**src/mcmc_inference/MCMC_temperature.py**

The project was originally designed to run either individual stellar cases or larger stellar samples using the functions provided by the modelling and inference modules. *The MCMC analysis produces posterior distributions and convergence plots for the inferred parameters.*


## Limitations

Several limitations identified during the original project remain relevant:

* Some stars cannot be processed because of missing or invalid catalogue data.
* MCMC convergence can be difficult or slow for some extreme stellar parameters.
* Very large stellar radii can be harder to recover reliably with the current initialization and sampling configuration.
* Low-temperature stars can approach the boundary of the available model grid.
* Extinction is currently treated using a simplified approach rather than being fully inferred for every analysis.

These limitations are useful considerations when applying the pipeline to larger or different datasets.

## Background

This project originated from my MSc dissertation:  

>"Measuring better stellar radius with Gaia: an application to transiting exoplanets"  
>  
>MSc Astronomy and Astrophysics  
>Faculty of Sciences, University of Porto (FCUP)

The original dissertation code has been reorganized and refactored here to make the data-processing and statistical-modelling workflow easier to understand, reproduce, and extend.
