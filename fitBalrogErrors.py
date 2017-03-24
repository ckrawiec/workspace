import numpy as np
from astropy.table import Table, join
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import glob

flux_col = 'AUTO'

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
    filename = '/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'.format(tab_num)
    
    print filename
    
    both = Table.read(filename)
    print "    # objects, joined: ", len(both)

    if i == 1:
        balrog = both
    else:
        balrog = join(balrog, both, join_type='outer')
    print "    # objects, total: ", len(balrog)
    
#trim to just galaxies
balrog = balrog[balrog['OBJTYPE']==1]
print "# objects, OBJTYPE==1: ", len(balrog)

def Gfunc(x, A, C, V):
    return A * np.exp(-0.5 * (x-C)*(x-C) / V)
    
def Efunc(x, A, B, C1, C2, V1, V2):
    return Gfunc(x, A, C1, V1) + Gfunc(x, B, C2, V2)


#error distribution bins
histbins = np.linspace(-20, 20, 10)

n_plot = 220
plt.subplots(2, 2, figsize=(15,10))

for b in 'GRIZ':
    flux0 = balrog['FLUX_0_'+b]
    size = balrog['HALFLIGHTRADIUS_0_'+b]
    
    flux = balrog['FLUX_'+flux_col+'_'+b]
    error = balrog['FLUXERR_'+flux_col+'_'+b]
    
    dig = np.digitize(flux0, flux_bins, right=True)
    edist = (flux-flux0)/error

    n_plot+=1
    for si in range(1, n_chunks+1):
        mask = np.where(dig==si)
        plt.subplot(n_plot)
        plt.hist(flux0[mask],
                 histtype='step', lw=2.,
                 bins=100, label=str(si))
        plt.xlabel('true flux')
        plt.xlim(0,5000)
    plt.title(b)
    plt.legend()
plt.show() 

f_digitized = np.digitize(balrog['FLUX_0_G'], flux_bins, right=True)
s_digitized = np.digitize(balrog['HALFLIGHTRADIUS_0_G'], size_bins, right=True)

for ibin in range(1, len(flux_bins)+1):
    for jbin in range(1, len(size_bins)+1):
        for band in 'GRIZ':
            error_func = (balrog['FLUX_'+flux_col+'_'+band]-balrog['FLUX_0_'+band]) / balrog['FLUXERR_'+flux_col+'_'+band]
            mask = ~np.isinf(error_func) & (np.abs(error_func) < 50) & (f_digitized==ibin) & (s_digitized==jbin) & (np.abs(balrog['FLUXERR_'+flux_col+'_'+band])<1000000.)
            y, bins, ps = plt.hist(error_func[mask], histtype='step', bins=30, normed=True)
            x = (bins[:-1]+bins[1:])/2
            plt.scatter(x, y)    
            #plt.title(str(ibin)+str(jbin))
        
            try:
                params, cov = curve_fit(Efunc, x, y)
                xplot = np.arange(x.min(),x.max(),0.1)
                Eargs = [xplot]+list(params)

                print ibin, jbin, band, params
                
                plt.plot(xplot, Efunc(*Eargs), lw=2., label=band)
            
                #plt.title(str(flux_bins[ibin]))

        
            except RuntimeError:
                continue
        plt.ylim(0, 0.4)
        plt.legend()
        plt.show()
    
