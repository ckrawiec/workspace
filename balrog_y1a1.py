import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join

truth_file = '/Users/Christina/DES/data/balrog_y1a1_truth_2m.fits' 
sim_file = '/Users/Christina/DES/data/balrog_y1a1_sim_2m_000001.fits'
nosim_file = '/Users/Christina/DES/data/balrog_y1a1_nosim_000001.fits'

#truth = fits.open(truth_file)[1].data
#sim = fits.open(sim_file)[1].data
nosim = fits.open(nosim_file)[1].data

truth = Table.read(truth_file)
sim = Table.read(sim_file)
joined = join(truth, sim)#, keys='balrog_index')

print len(truth), len(sim), len(joined)
plt.scatter(joined['FLUX_0_I'], joined['FLUX_DETMODEL_I'], edgecolor='none', c=joined['Z'])
plt.colorbar(label='z')
plt.xlabel('FLUX_0_I')
plt.ylabel('FLUX_DETMODEL_I')
plt.ylim(-1000, 200000)
plt.xlim(0, 200000)
plt.savefig('balrog_test')
plt.close()

plt.hist(truth['Z'], bins=1000)
plt.xlabel('truth z')
plt.savefig('balrog_truth_hist_z')
plt.close()

#look at distributions
#compare to DES cosmos
#MAKE TABLE, JUST MAKE ONE
#DES y1a1 cosmos - is there already a list of coadd ids for cosmos? (maybe only sv)
#if not, can just match by ra/dec

#in one of balrog tables, does it tell you if it was not detected, or is it just not included
