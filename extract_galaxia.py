import numpy as np
import glob
import time
from astropy.io import fits
from astropy.table import Table

n_files = -1
home_dir = '/home/ckrawiec/'

#galaxia files
gals = glob.glob(home_dir+'DES/data/Nside128.11*.fits')

output = home_dir+'DES/data/DES_NSIDE128_11.fits'#.format(n_files)

g,r,i,z = [],[],[],[]
glon, glat = [],[]

start = time.time()
for gal in gals[:n_files]:
    gal_tab = fits.open(gal, memmap=False)[1].data
    
    g.append(gal_tab['sdss_g'])
    r.append(gal_tab['sdss_r'])
    i.append(gal_tab['sdss_i'])
    z.append(gal_tab['sdss_z'])

    glon.append(gal_tab['glon'])
    glat.append(gal_tab['glat'])
    
end = time.time()

#make new table
tab = Table()

tab['glon'] = np.hstack(glon)
tab['glat'] = np.hstack(glat)

tab['sdss_g'] = np.hstack(g)
tab['sdss_r'] = np.hstack(r)
tab['sdss_i'] = np.hstack(i)
tab['sdss_z'] = np.hstack(z)

tab.write(output)
print "completed in {}s".format(end-start)
