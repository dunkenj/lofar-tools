import os
import numpy as np
from scipy.special import erfinv
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits

from astropy.coordinates import SkyCoord
from astropy.stats import median_absolute_deviation, sigma_clipped_stats
import astropy.visualization as vis 
from astroquery.ibe import IbeClass

from ccdproc import CCDData, wcs_project

import urllib
import urlparse


def load_image(filename):
    try:
        data = CCDData.read(filename)
    except(ValueError):
        try:
            data = CCDData.read(filename, unit = u.dyn)
        except:
            raise
    return data

def split_url(url):
    """
    Extract filename from full url of cutout image
    
    Paramters
    ---------
    url : string
        url from which to extract filename
        
    Returns
    -------
    filename : string
        extracted filename
    
    """
    url_split = urlparse.urlsplit(url)
    head, filename = os.path.split(url_split.path)
    return filename

def dl_wise(coordinates, ang_size, band=1, clobber=True):
    """
    Get cutouts of ALLWISE images
        
    
    Parameters
    ----------
    coordinates : astropy.coordinates object
        Central coordinates of desired cutout
        
    ang_size : astropy.quantity object
        Angular diamater of desired cutout (max = 5 degrees)
    
    band : integer, default = 1
        WISE band to get cutout for. WISE filters correspond to rest-frame
        central wavelengths of W1 = 3.4um, W2 = 4.6um, W3 = 12.0um and 
        W4 = 22.0um.
        
    clobber : bool, default = True
        Overwrite existing filenames when generating cutouts
        
    Returns
    -------
    
        
    """
    
    pixel_scale= 1.3333333333*u.arcsec
    size_pix = int(ang_size/pixel_scale)
    
    Ibe = IbeClass()
    out = Ibe.query_region_sia(coordinates, mission='wise', 
                               dataset='allwise', 
                               table='p3am_cdd', width = ang_size,
                               intersect='CENTER',
                               most_centered=True)
                           
    query_string ='?center={0},{1}'\
                  'deg&size={2}arcsec'.format(coordinates.ra.to(u.deg).value,
                                           coordinates.dec.to(u.deg).value,
                                           ang_size.to(u.arcsec).value)

    obj_url = out['sia_url'][band-1]
    filename = split_url(obj_url)
    if os.path.isfile(filename) and clobber:
        os.remove(filename)
        
    test = urllib.urlretrieve(obj_url+query_string,
                              filename=filename)
    print(obj_url+query_string)

    return filename

def get_wise(coordinates, ang_size, band=1, clobber=True):
    """
    Get cutouts of ALLWISE images
        
    
    Parameters
    ----------
    coordinates : astropy.coordinates object
        Central coordinates of desired cutout
        
    ang_size : astropy.quantity object
        Angular diamater of desired cutout (max = 5 degrees)
    
    band : integer, default = 1
        WISE band to get cutout for. WISE filters correspond to rest-frame
        central wavelengths of W1 = 3.4um, W2 = 4.6um, W3 = 12.0um and 
        W4 = 22.0um.
        
    clobber : bool, default = True
        Overwrite existing filenames when generating cutouts
        
    Returns
    -------
    
        
    """
    filename = dl_wise(coordinates, ang_size, band=1, clobber=True)
    
    image = load_image(filename)
        
    mad = median_absolute_deviation(image)
    sigma = mad#/np.sqrt(2)*erfinv(2*0.75-1)

    image = image.subtract(np.median(image)*image.unit)
    image = image.divide(sigma)

    return image #, out

def get_decals(coordinates, size, band='g', 
               pix_scale=0.262*u.arcsec, clobber=True):
               
    size_pix = np.round((size / pix_scale),0)+1
    if size_pix > 512:
        raise Exception()
    
    base_url='http://legacysurvey.org/viewer/fits-cutout-decals-dr2?'
    query_string = 'ra={0}&dec={1}&pixscale={2}'\
                   '&size={3:.0f}&bands={4}'.format(
                                              coordinates.ra.to(u.deg).value,
                                              coordinates.dec.to(u.deg).value,
                                              pix_scale.to(u.arcsec).value,
                                              size_pix,
                                              band)

    filename = 'test.fits'
    if os.path.isfile(filename) and clobber:
        os.remove(filename)
    print base_url+query_string
    test = urllib.urlretrieve(base_url+query_string,
                              filename=filename)

    image = load_image(filename)        
    mad = median_absolute_deviation(image)
    sigma = mad#/np.sqrt(2)*erfinv(2*0.75-1)

    image = image.subtract(np.median(image)*image.unit)
    image = image.divide(sigma)


    return image #, out

def scale_rgb(red, green, blue, stretch_coeff,
             pmin=0, pmax=98.):
    """ Perform interval and stretch for RGB
    
    Parameters
    ----------
    red : CCDData object
        Red filter image data, parsed with load_image
    green : CCDData object
        Green filter image data, parsed with load_image    
    blue : CCDData object
        Blue filter image data, parsed with load_image
    stretch_coeff : tuple, (1,3)
        Power-law stretch coefficients for (red, green, blue)
        filters respectively.
    pmin : float, default = 0
        The lower percentile below which to ignore pixels. Passed to 
        AsymmetricPercentileInterval function.
    pmax : float, default = 98
        The upper percentile above which to ignore pixels. Passed to 
        AsymmetricPercentileInterval function.    
        
    Returns
    -------
    rgb : np.array, N x N x 3
    
    """
    r_stretch = vis.stretch.PowerStretch(stretch_coeff[0])
    g_stretch = vis.stretch.PowerStretch(stretch_coeff[1])
    b_stretch = vis.stretch.PowerStretch(stretch_coeff[2])
    
    interval = vis.AsymmetricPercentileInterval(pmin, pmax)
    
    rgb = np.array([r_stretch(interval(red)),
                    g_stretch(interval(green)),
                    b_stretch(interval(blue))]).T
    return rgb