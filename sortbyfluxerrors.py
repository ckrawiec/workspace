import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table

t_file ='/home/ckrawiec/DES/data/y1a1_gold_flux_detmodel_MC1.fits'
#'/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel.fits'

filters = ['G','R','I','Z']
N_split = 100

data = Table.read(t_file)

data_name = os.path.splitext(t_file)[0]

fluxerr = (data['FLUXERR_DETMODEL_'+filter] for filter in filters)

errs = np.array(zip(*fluxerr))
sigmas = np.sqrt(np.sum(errs**2., axis=1))
ids = sigmas.argsort()

idsplit = np.array_split(ids, N_split)

print "Length of each split:"
s=0
for i in idsplit:
    print len(i)
    s+= len(i)

print "Total objects:", s
print "Total in full table: ", len(data)

tabs=[data[id] for id in idsplit]

del data

tabnames = [data_name+'_'+''.join(filters)+'_fluxerrgrp{0:0>3}.fits'.format(ir+1) for ir in range(N_split)]
for t,n in zip(tabs, tabnames):
    if os.path.exists(n):
        print "File {} already exists, delete before running again".format(n)
    else:
        t.write(n)

for filter in filters:
    ymins = np.array([np.percentile(tab['FLUXERR_DETMODEL_'+filter], 10) for tab in tabs])
    ymaxs = np.array([np.percentile(tab['FLUXERR_DETMODEL_'+filter], 90) for tab in tabs])
    means = np.array([tab['FLUXERR_DETMODEL_'+filter].mean() for tab in tabs])
    plt.errorbar(range(N_split), means, 
                 yerr=[means-ymins,ymaxs-means], fmt='o-', label=filter)
plt.ylim(1,1100000)
plt.xlabel('flux error group')
plt.ylabel('mean fluxerr_detmodel')
plt.title('errorbars = 10th & 90th percentile')
plt.yscale('log')
plt.legend(loc='best')
plt.savefig('{}_{}_fluxerrgrp_means'.format(data_name,''.join(filters)))
plt.close()

for filter in filters:
    ymins = np.array([np.percentile(tab['FLUXERR_DETMODEL_'+filter],10) for tab in tabs])
    ymaxs = np.array([np.percentile(tab['FLUXERR_DETMODEL_'+filter],90) for tab in tabs])
    medians = np.array([np.median(tab['FLUXERR_DETMODEL_'+filter]) for tab in tabs])
    plt.errorbar(range(N_split), medians, 
                 yerr=[medians-ymins,ymaxs-medians], fmt='o-', label=filter)
plt.ylim(10,1100)
plt.xlabel('flux error group')
plt.ylabel('median fluxerr_detmodel')
plt.yscale('log')
plt.title('errorbars = 10th & 90th percentile')
plt.legend(loc='best')
plt.savefig('{}_{}_fluxerrgrp_medians'.format(data_name,''.join(filters)))
plt.close()

