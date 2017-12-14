# Questions to answer:
# - how close do i need them to be?
# - check mag differences? see how close compared to surrouning x objects?
# - use fraction of closest objects that were detected vs. magnified ones?
# - is there another way to match them?
import matplotlib.pyplot as plt
import numpy as np
import pandas
import sys
import time
from astropy.table import Table, join
from astropy.io import fits
from scipy.spatial import ckdtree
from scipy.interpolate import griddata
from matplotlib import cm

#frequency of objects to query tree with for now
freq = 100

# Read Data
#Balrog truth table - eventually all tables
obrog = fits.open('/Users/Christina/DES/data/balrog_sva1_tab01_TRUTH_fluxes.fits')[1].data

#DFULL in COSMOS region
odfull = fits.open('/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits')[1].data
#all DFULL
adfull = fits.open('/Users/Christina/DES/data/y1a1_gold_dfull_radii.fits')[1].data

#Balrog SV sim catalogs - flux radii and auto values
rsim = fits.open('/Users/Christina/DES/data/balrog_sva1_tab01_SIM_flux_radii.fits')[1].data
#sim = fits.open('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab01_SIM.fits')[1].data

#cut out low fluxes - remove those below ~noise level
cut = 40
flux_mask = (odfull['FLUX_AUTO_G'] > cut) & (odfull['FLUX_AUTO_R'] > cut) & (odfull['FLUX_AUTO_I'] > cut) & (odfull['FLUX_AUTO_Z'] > cut)

#choose just galaxies - could change modest_class to something else
brog_gals = obrog['OBJTYPE']==1
dfull_gals = odfull['MODEST_CLASS']==1
brog = obrog[brog_gals]
dfull = odfull[dfull_gals & flux_mask]

#join truth galaxies and sim radii
bsim = join(brog, rsim)
#join all objects and sim radii
obsim = join(obrog, rsim)

#other masks
ostars = odfull['MODEST_CLASS']==0
astars = adfull['MODEST_CLASS']==0
agals = adfull['MODEST_CLASS']==1
bstars = obsim['OBJTYPE']==3 #or MODEST_CLASS?

# Check Size Scaling
# in arcsec
dfull_size = dfull['FLUX_RADIUS_I']
brog_size = brog['HALFLIGHTRADIUS_0']

# Relating Balrog & DFULL Sizes
# Use stars to get average PSF size around DFULL galaxies 

#make tree of dfull stars
sys.stderr.write("stars in all dfull: {}\n".format(len(adfull[astars])))
sys.stderr.write("galaxies in all dfull: {}\n".format(len(adfull[agals])))

sys.stderr.write('creating tree...\n')
star_tree = ckdtree.cKDTree(zip(adfull[astars]['RA'],
                        adfull[astars]['DEC']))

gal_pos = zip(adfull[agals]['RA'][::freq], adfull[agals]['DEC'][::freq])
#360 arcsec
close = star_tree.query_ball_point(gal_pos, r=0.1)
sys.stderr.write("number of matches for first object: {}\n".format(len(close[0])))

#make tree of balrog stars
sys.stderr.write("stars in sim balrog: {}\n".format(len(obsim[bstars])))
sys.stderr.write("galaxies in sim balrog: {}\n".format(len(bsim)))
sys.stderr.write('creating tree...\n')

brog_star_tree = ckdtree.cKDTree(zip(obsim['RA'][bstars], obsim['DEC'][bstars]))
sim_pos = zip(bsim['ALPHAMODEL_J2000_G'][::freq], bsim['DELTAMODEL_J2000_G'][::freq])
brog_close = brog_star_tree.query_ball_point(sim_pos, r=0.1)

#calculate average star size around each galaxy
sys.stderr.write('calculating average star radii...\n')

start = time.time()
avg_psf = np.array([np.median(adfull['FLUX_RADIUS_I'][astars][close[i]]) for i in range(len(close))])
brog_avg_psf = np.array([np.median(obsim['FLUX_RADIUS_I'][bstars][brog_close[i]]) for i in range(len(brog_close))])
end = time.time()

sys.stderr.write('time to calculate radii = {}s\n'.format(end-start))

df = pandas.DataFrame({'brog_psf': brog_avg_psf, 
                      'brog_rad': bsim['FLUX_RADIUS_I'][::freq],
                      'brog_hlr': bsim['HALFLIGHTRADIUS_0'][::freq]})

interp = griddata(zip(df['brog_rad'], df['brog_psf']),
                  df['brog_hlr'],
                  zip(adfull['FLUX_RADIUS_I'][astars], avg_psf))

plt.hist(df['brog_hlr'], histtype='step', normed=True, bins=100,
         label='brog hlr')
plt.hist(interp[~np.isnan(interp)], histtype='step', normed=True, bins=100,
         label='interp')
plt.legend(loc='best')
plt.show()

tab = Table()
tab['flux_radius_i'] = df['brog_rad']
tab['avg_flux_radius_i'] = df['brog_psf']
tab['hlr'] = df['brog_hlr']
sys.stderr.write('writing to table...\n')
tab.write('/Users/Christina/DES/data/balrog/sva1/balrog_tab01_avg_star_fluxradiusi_0.1deg.fits')
