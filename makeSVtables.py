import numpy as np
from astropy.table import Table, Column
import healpy as hp

home_dir = '/Users/Christina/'

sva1_gold = Table.read(home_dir+'DES/data/sva1_gold_detmodel_gals.fits')
cosmos = Table.read(home_dir+'COSMOS/data/COSMOS2015_Laigle+_v1.1.fits')

print "loaded files"

region_file = home_dir+'DES/data/sva1_gold_r1.0_goodregions_04_n4096.fits.gz'

hpmap = hp.read_map(region_file, nest=True)
nside = hp.npix2nside(hpmap.size)

theta = (90.0 - sva1_gold['dec'])*np.pi/180.
phi = sva1_gold['ra']*np.pi/180.
pix = hp.ang2pix(nside, theta, phi, nest=True)
good, = np.where(hpmap[pix] == 1)

print "good regions found"

h = esutil.htm.HTM(10)
h.match(sva1_gold['ra'][good,], sva1_gold['dec'][good,],
        cosmos['alpha_j2000'], cosmos['delta_j2000'],
        radius=1./3600,
        file=home_dir+'DES/data/match_sva1_gold_cosmos_gals_1arcsec')
m = h.read(home_dir+'DES/data/match_sva1_gold_cosmos_gals_1arcsec')

gold_m, cosmos_m, merr = zip(*m)

print "matched"

no_cosmos = list(set(range(0,len(sva1_gold[good,])))-set(gold_m))

new_sv = sva1_gold[good,][no_cosmos]

print "made new SV table"

new_sv.write(home_dir+'DES/data/sva1_gold_detmodel_MC1_good_regions_no_cosmos.fits')
del new_sv

print "wrote new SV table"

new_cosmos = sva1_gold[good,][gold_m] 
del sva1_gold

print "deleted old SV table"

new_cosmos.add_column(Column(name='photoz', data=cosmos['PHOTOZ'][cosmos_m]))
new_cosmos.add_column(Column(name='match_err', data=merr))

print "made new cosmos table"

new_cosmos.write(home_dir+'DES/data/sva1_gold_detmodel_MC1_good_regions_cosmos.fits')
