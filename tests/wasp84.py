"""
WASP-84 example
----------------

Runs the full stellar-radius inference pipeline for WASP-84
using the MCMC temperature model.
"""

from src.mcmc_inference.MCMC_temperature import star_tester_main
from src.data_retrieval.query_sweetcat import query_sweetcat
from astropy.constants import R_sun

# Query SWEET-cat database for stellar parameters
star_name = "WASP-84"
sweetcat_data = query_sweetcat(star_name)
star = sweetcat_data[0]

# Extract stellar parameters
Teff = star["Teff"]
Teff_unc = star["e_Teff"]

log_g = star["logg"]
log_unc = star["e_logg"]

metallicity = star["__Fe_H_"]
met_unc = star["e__Fe_H_"]

# Reference radius used for comparison
reference_radius = star["Radius-t"] * R_sun
reference_radius_unc = star["e_Radius-t"] * R_sun

Ebv = 0.020


# MCMC settings
nwalkers = 32
npoints = 1000


# Run the full analysis
radius, radius_uncertainty, sampler = star_tester_main(
    star_name,
    Teff,
    Teff_unc,
    log_g,
    log_unc,
    metallicity,
    met_unc,
    Ebv,
    reference_radius,
    nwalkers,
    npoints,
)


# Print final result
print("\n--- WASP-84 RESULT ---")
print(f"Inferred radius: {radius.value:.3f} R_sun")
print(f"Uncertainty:     ±{radius_uncertainty:.3f} R_sun")

reference_radius_value = reference_radius.to(R_sun).value
print(f"Reference radius: {reference_radius_value:.3f} R_sun")

difference = (abs(radius.value - reference_radius_value)  / reference_radius_value * 100)

print(f"Difference:       {difference:.2f}%")
