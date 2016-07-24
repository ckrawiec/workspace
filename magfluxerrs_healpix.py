import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import pandas
from astropy.io import fits

map_file = "/Users/Christina/DES/sva1_gold_1.0.2-4_nside4096_nest_i_detmodel_weights.fits.gz"
sva1_gold = fits.open('/Users/Christina/DES/data/sva1_gold_detmodel_gals.fits')[1].data

ra, dec = sva1_gold['ra'], sva1_gold['dec']

theta = (90.0 - dec)*np.pi/180.
phi = ra*np.pi/180.

magnside = 512

magpix = hp.ang2pix(magnside, theta, phi, nest=True)

maghpmap = np.zeros(hp.nside2npix(magnside), dtype=np.double)

#maghpmap[magpix] = sva1_gold['mag_detmodel_i']

print magpix.min(), magpix.max(),len(maghpmap),len(sva1_gold), len(theta), len(magpix)

df = pandas.DataFrame( {'mag': sva1_gold['mag_detmodel_i'],
                        'magerr': sva1_gold['magerr_detmodel_i'],
                        'pixel': magpix} )

gpix = df.groupby(['pixel'])
plt.scatter(gpix['mag'].mean(),gpix['magerr'].mean())
plt.show()

#hpmap = hp.read_map(map_file, nest=True)
#nside = hp.npix2nside(hpmap.size)
#pix = hp.ang2pix(nside,theta,phi,nest=True)
#good, = np.where(hpmap[pix] > 0)

#hp.mollview(hpmap, nest=True)

#N = 10000000
#plt.scatter(ra[good][:N], dec[good][:N], c=hpmap[pix][good][:N], vmin=23.4, vmax=23.6, edgecolor='none', s=4.)
#plt.colorbar(label='limiting magnitude, detmodel i')
#plt.xlim(104,106)
#plt.ylim(-56.5,-56)

#hpmap.size

#len(hpmap[np.where(hpmap>0)])/4096

#ang = hp.pix2ang(4096, pix)

#mynside = 1024
#mypix = hp.ang2pix(mynside, theta, phi, nest=True)
#healpix_map = np.zeros(hp.nside2npix(mynside), dtype=np.double)
#healpix_map[mypix] = ra

#hp.mollview(healpix_map, nest=True)

#magerrhpmap = np.zeros(hp.nside2npix(magnside), dtype=np.double)
#magerrhpmap[magpix] = sva1_gold['magerr_detmodel_i']
#fluxhpmap = np.zeros(hp.nside2npix(magnside), dtype=np.double)
#fluxerrhpmap = np.zeros(hp.nside2npix(magnside), dtype=np.double)
#fluxhpmap[magpix] = sva1_gold['flux_detmodel_i']
#fluxerrhpmap[magpix] = sva1_gold['fluxerr_detmodel_i']

#Nmag=10000
#plt.scatter(maghpmap[magpix][:Nmag], magerrhpmap[magpix][:Nmag], edgecolor='none', s=4.)
#plt.xlim(18,30)
#plt.ylim(0.,1.)
#plt.xlabel('mag_detmodel_i')
#plt.ylabel('magerr_detmodel_i')
#plt.title('NSIDE='+str(magnside))



"""
plt.scatter(fluxhpmap[magpix][:Nmag], fluxerrhpmap[magpix][:Nmag], edgecolor='none', s=4.)
plt.xlim(18,30)
plt.ylim(0.,1.)
plt.xlabel('flux_detmodel_i')
plt.ylabel('fluxerr_detmodel_i')
plt.title('NSIDE='+magnside)
"""


