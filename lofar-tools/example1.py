import numpy as np
import astropy.units as u


asize =  60*u.arcsec
test_pos = SkyCoord(244.6461*u.deg, 7.4106*u.deg)
wise = get_wise(test_pos, asize)
pscale = (wise.header['PXSCAL1']*u.arcsec)
pscale = 0.5*u.arcsec


g = get_decals(test_pos, asize, band='g',
                      pix_scale=pscale)
r = get_decals(test_pos, asize, band='r',
                      pix_scale=pscale)
z = get_decals(test_pos, asize, band='z',
                      pix_scale=pscale)

wise2 = wcs_project(wise, g.wcs, g.shape)

imstretch = vis.stretch.LogStretch()

rgb = scale_rgb(z, r, g, np.array([1, 1, 1])*1., pmin=1, pmax=98)

rgb_wise = scale_rgb(wise2, z, g, np.array([0.66, 1, 1])*1., pmin=1, pmax=98)

Fig, Ax = plt.subplots(1,2)
Ax[0].imshow(rgb, interpolation='nearest')
Ax[1].imshow(rgb_wise, interpolation='nearest')
plt.show()