import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from myutils import fitsstack, match

tables = ['/Users/Christina/DES/data/balrog/y1a1/balrog_y1a1_truth_objtype1_fluxg_000001.fits',
          '/Users/Christina/DES/data/balrog/y1a1/balrog_y1a1_truth_objtype1_fluxg_000002.fits',
          '/Users/Christina/DES/data/balrog/y1a1/balrog_y1a1_truth_objtype1_fluxg_000003.fits',
          '/Users/Christina/DES/data/balrog/y1a1/balrog_y1a1_truth_objtype1_fluxg_000004.fits']

hdu = fitsstack(tables)
data = hdu.data

mus = np.arange(1.01, 1.05, 0.01)

def plothist(col):
    plt.hist(data[col], histtype='step', bins=1000)
    plt.xlabel(col)
    plt.savefig('/Users/Christina/DES/data/balrog/y1a1/'+col+'_hist')
    plt.close()

plothist('flux_0_g')
plothist('flux_noised_g')
plothist('flux_noiseless_g')

flux = data['flux_0_g']

print "Entries in table: " + str(len(flux))

flux.sort()
plt.hist(flux[:10000000], bins=1000, label='flux', histtype='step')

for mu in mus:
    print "working on mu={}...".format(mu)
    newflux = flux*mu
    f, nf = match(flux, newflux)
    print "    {} matches found".format(len(f))
    plt.hist(newflux[:10000000], bins=1000, label="mu={}".format(mu), histtype='step')

plt.legend(loc='best')
plt.savefig('/Users/Christina/DES/data/balrog/y1a1/flux_mu_hist')
plt.close()

print flux[f][0]
print newflux[nf][0]
