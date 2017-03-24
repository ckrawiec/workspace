import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.table import Table, join

#TODO: look at objects where flux - flux0 ~=0, decide how to remove
#      figure out how to match gaussian without checking first
#      check difference between flux0/flux_noiseless for these plots
#      size columns
#      reasonable bin sizes/edges
#      look at medians in different percentiles
#      fit a gaussian/other distribution?

flux_cols = ['AUTO', 'DETMODEL', 'APER_5', 'APER_10']
size_col = 'HALFLIGHTRADIUS_0_{}'
true_flux_col = 'NOISELESS'

bands = 'G'

t = Table.read('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_tab01_TRUTH.fits')

#slices in which to plot
flux_bins = np.logspace(0, 5, 4)
print flux_bins
size_bins = np.linspace(0, 2, 4)
print size_bins

#subplots
start_plotn = 220


for flux_col in flux_cols:
    a = Table.read('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_{}_tab01_SIM.fits'.format(flux_col))

    at = join(t, a)

    f, axs = plt.subplots(3, 3, figsize=(10,10))
    plt.xlabel('(FLUX_{} - FLUX_{}) / FLUXERR_{}'.format(flux_col, true_flux_col, flux_col))
    
    for b in bands:
        plotn = start_plotn
        size = at[size_col.format(b)]
        flux = at['FLUX_{}_{}'.format(flux_col, b)]
        flux_err = at['FLUXERR_{}_{}'.format(flux_col, b)]
        true_flux = at['FLUX_{}_{}'.format(true_flux_col, b)]

        flux_digitized = np.digitize(true_flux, flux_bins)
        size_digitized = np.digitize(size, size_bins)

        for i in range(len(flux_bins)):
            plotn +=1

            ax = plt.subplot(plotn)
            for j in range(len(size_bins)):
                mask = (flux_digitized==i) & (size_digitized==j)
                nerr_dist = (flux[mask] - true_flux[mask]) / flux_err[mask]
                nerr_dist = nerr_dist[~np.isnan(nerr_dist) & ~np.isinf(nerr_dist)]
                nerr_dist = nerr_dist[(np.abs(nerr_dist) < 20) & (np.abs(nerr_dist) > 0.02)]
                h = ax.hist(nerr_dist,
                            bins=50, histtype='step',
                            label=str(size_bins[j]),
                            normed=True)
            #ax.legend()
            ax.grid()
            xticklabels = ax.get_xticklabels()
            yticklabels = ax.get_yticklabels()
            ax.set_xticklabels(xticklabels, fontsize=10)
            ax.set_yticklabels(yticklabels, fontsize=10)
            x = np.arange(-20,20,40./50.)
            y = norm.pdf(x)
            #y = np.exp(-0.5*x*x) #/ np.sqrt(2*np.pi) * (100./100000)
            plt.plot(x,y,'k--',lw=2)
            #ax.set_title('flux in {}'.format([np.min(flux[flux_digitized==i]),np.max(flux[flux_digitized==i])]))
        #topbins = np.argsort(h[1][:-1])[-100:]
        #ctr = np.median(h[0][topbins])
        #print ctr ###
    #y = len(at)* 500* np.exp(-0.5*x*x)
    plt.subplots_adjust(hspace=0.23, wspace=0.23)
    plt.xlim(-10,10)
    plt.ylim(0, 0.4)
    #plt.legend()
    #plt.grid()
    plt.savefig('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_{}_flux{}err_fluxbins_dist'.format(flux_col, true_flux_col))
    plt.close()
