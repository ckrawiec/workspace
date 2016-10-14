import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join

balrog = Table.read('/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_size.fits')

mu = 1.02

mag = Table.read('/home/ckrawiec/DES/magnification/balrog_mag_02.fits')

filters = ['G','R','I','Z','Y']

new = join(balrog, mag)

del balrog
del mag

"""
plt.hist(new['D_MAG1.01'], bins=np.logspace(0,5,100))
plt.xlabel('distance to nearest neighbor in truth tree')
plt.xscale('log')
plt.title('$\mu=1.01$')
plt.savefig('/home/ckrawiec/DES/magnification/hist_balrog_d_mag_01')
plt.close()

match = np.where(new['INDEX_MAG1.01']==new['BALROG_INDEX'])
gal = np.where(new['OBJTYPE']==1)

print # where nearest neighbor is same object: {}.format(len(match[0]))

print # where objtype=1 (galaxy): {}.format(len(gal[0]))
"""

for filter in filters:

    newflux = new['FLUX_NOISELESS_'+filter]*mu
    foundflux = new['FLUX_NOISELESS_'+filter+'_MAG'+str(mu)]
    diff = foundflux - newflux

    plt.hist(diff, bins=5000, histtype='step', label=filter)

plt.legend(loc='best')
plt.xlabel('$found flux - \mu * flux$')
plt.title('$\mu=$'+str(mu))
#plt.savefig('/home/ckrawiec/DES/magnification/balrog_diff_hist_lin_mag_'+str(mu)[2:])
plt.show()

newsize = new['HALFLIGHTRADIUS_0']
foundsize = new['HALFLIGHTRADIUS_0_'+str(mu)]
diffsize = foundsize - newsize
plt.hist(diffsize, bins=5000, histtype='step')
plt.xlabel('$found hlr - \mu * hlr$')
#plt.savefig('/home/ckrawiec/DES/magnification/balrog_diff_hlr_hist_mag'+str(mu)[2:])
plt.show()
