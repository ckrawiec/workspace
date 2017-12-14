import numpy as np
import glob
from astropy.io import fits
from astropy.table import Table

#galaxia files
gals = glob.glob('/Users/Christina/DES/data/DES_NSIDE128/Nside128.*.fits')

output = '/Users/Christina/DES/data/DES_NSIDE128_all.fits'

g,r,i,z = [],[],[],[]
glon, glat = [],[]

for gal in gals:
    gal_tab = fits.open(gal, memmap=False)[1].data
    
    g.append(gal_tab['sdss_g'])
    r.append(gal_tab['sdss_r'])
    i.append(gal_tab['sdss_i'])
    z.append(gal_tab['sdss_z'])

    glon.append(gal_tab['glon'])
    glat.append(gal_tab['glat'])
    
#make new table
tab = Table()

tab['glon'] = np.hstack(glon)
tab['glat'] = np.hstack(glat)

tab['sdss_g'] = np.hstack(g)
tab['sdss_r'] = np.hstack(r)
tab['sdss_i'] = np.hstack(i)
tab['sdss_z'] = np.hstack(z)

tab.write(output)
