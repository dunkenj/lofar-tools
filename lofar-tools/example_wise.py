import numpy as np
import astropy.units as u
from cutouts import *
from astropy.coordinates import SkyCoord


asize =  60*u.arcsec # Image size
sky_positions = SkyCoord([[244.2, 7.4],
                          [245.6, 7.42]], 
                          unit=[u.deg, u.deg])

for i, coordinates in enumerate(sky_positions):
    # w1_path = filename of downloaded image
    w1_path = dl_wise(coordinates, asize, band=1) # 3.4um
    w2_path = dl_wise(coordinates, asize, band=2) # 4.6um

    
    # Rename files to some ID number convenience (if needed)
    idnumber = i+1 # or whatever radio ID etc. you want
    os.rename(w1_path, '{0}_WISE_W1.fits'.format(idnumber))
    os.rename(w2_path, '{0}_WISE_W2.fits'.format(idnumber))

