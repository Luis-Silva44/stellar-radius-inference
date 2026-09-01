import numpy as np
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from uncertainties import ufloat
import math

from src.data_retrieval.auxiliary_functions import retrieve_gaia_id, mag_to_flux, vizier_coords

def wise_values(star_name):
    ''' Function to retrieve W1 and W2 magnitudes from WISE survey and turn them into fluxes. Also checks for availability of each band'''
    
    Vizier.ROW_LIMIT = -1
    _, star_name = retrieve_gaia_id(star_name)
    wise_catalog = 'II/311/wise'
    ra, dec = vizier_coords(star_name)
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    wise_data = Vizier.query_region(coords, radius=10* u.arcsec , catalog=wise_catalog)

    if wise_data:
        W1_mag = ufloat(wise_data[0][0]['W1mag'], wise_data[0][0]['e_W1mag'])
        W2_mag = ufloat(wise_data[0][0]['W2mag'], wise_data[0][0]['e_W2mag'])
        
        W1_flux = mag_to_flux(W1_mag,'W1')
        W2_flux= mag_to_flux(W2_mag,'W2')

        unit = 1 * u.watt / u.um / u.cm**2
        wise_flux = np.array([W1_flux, W2_flux]) * unit
        wise_check = ['W1','W2']
        
    else:
        unit = 1 * u.watt / u.um / u.cm**2
        wise_check = []
        wise_flux = [] * unit
        print('No WISE data found')
    
    return wise_flux, wise_check
