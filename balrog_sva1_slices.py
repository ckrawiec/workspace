import numpy as np
from astropy.table import Table, join
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import glob

flux_col = 'AUTO'

tab_dir = '/Users/Christina/DES/data/balrog/sva1/'

#flux bins
n_chunks = 2
flux_bins = np.logspace(2, 3, n_chunks)

#size bins
n_chunks = 2
size_bins = np.logspace(-1, 0, n_chunks)

for i in range(1,3):
    if i == 9:
        continue
    tab_num = str(i).zfill(2)
    filename = tab_dir+'balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'.format(tab_num)
    print filename
    balrog = Table.read(filename)
    print "    # objects: ", len(balrog)

    for ibin in range(1, len(flux_bins)+1):
        for jbin in range(1, len(size_bins)+1):
            balrog_new = Table(balrog)
            for band in 'GRIZ':
                s_digitized = np.digitize(balrog_new['HALFLIGHTRADIUS_0_'+band], size_bins, right=True)
                f_digitized = np.digitize(balrog_new['FLUX_0_'+band], flux_bins, right=True)
                error_func = (balrog_new['FLUX_'+flux_col+'_'+band]-balrog_new['FLUX_0_'+band]) / balrog_new['FLUXERR_'+flux_col+'_'+band]
                mask = ~np.isinf(error_func) & (np.abs(error_func) < 50) & (f_digitized==ibin) & (s_digitized==jbin) & (np.abs(balrog_new['FLUXERR_'+flux_col+'_'+band])<1000000.) & (balrog_new['FLUX_'+flux_col+'_'+band]<5000)
                balrog_new = balrog_new[mask]
            print "    # object in flux bin {}, size bin {}: {}".format(ibin, jbin, len(balrog_new))
            balrog_new.write(tab_dir+'balrog_sva1_{}_tab{}_SIM_TRUTH_zp_corr_fluxes_fluxbin{}_sizebin{}.fits'.format(flux_col, tab_num, ibin, jbin))
            
