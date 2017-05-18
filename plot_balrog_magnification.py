import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join

mu = 1.04
mag_file = '/home/ckrawiec/DES/magnification/balrog_y1a1_truth_noiseless_flux_size_objtype1_80M_mag_04.fits'
filters = ['G','R','I','Z','Y']

def makeplots():
    fluxdiff()
    sizediff()
    Nvsdiffthresh()
    
balrog_file = '/home/ckrawiec/DES/data/balrog_y1a1_truth_noiseless_flux_size_objtype1.fits'
sim_file = '/home/ckrawiec/DES/data/balrog_y1a1_sim_truth_flux_detmodel.fits'

mag = Table.read(mag_file)
balrog = Table.read(balrog_file)
new = join(balrog, mag)

del balrog
del mag

def fluxdiff():
    for filter in filters:

        newflux = new['FLUX_NOISELESS_'+filter]*mu
        foundflux = new['FLUX_NOISELESS_'+filter+'_MAG'+str(mu)]
        diff = foundflux - newflux
        
        plt.hist(diff, bins=5000, histtype='step', label=filter)

        plt.legend(loc='best')
        plt.xlabel('$found flux - \mu * flux$')
        plt.title('$\mu=$'+str(mu))

    plt.xlim(-200,200)
    plt.savefig('/home/ckrawiec/DES/magnification/balrog_diff_hist_mag_'+str(mu)[2:])
    plt.show()

def sizediff():
    newsize = new['HALFLIGHTRADIUS_0'] * np.sqrt(mu) 
    foundsize = new['HALFLIGHTRADIUS_0_MAG'+str(mu)]
    diffsize = foundsize - newsize
    plt.hist(diffsize, bins=5000, histtype='step')
    plt.xlabel('$found hlr - \mu^{1/2} * hlr$')
    plt.xlim(-0.02,0.01)
    plt.savefig('/home/ckrawiec/DES/magnification/balrog_diff_hlr_hist5000_mag_'+str(mu)[2:])
    plt.close()

def Nvsdiffthresh():
    thresholds = np.arange(0,10,0.1)

    for filter in filters:
        newflux = new['FLUX_NOISELESS_'+filter] * mu
        foundflux = new['FLUX_NOISELESS_'+filter+'_MAG'+str(mu)]
        diff = foundflux - newflux

        N = [len(np.where(np.abs(diff) < thresh)[0]) for thresh in thresholds]
        plt.plot(thresholds, N, label=filter)

    plt.xlabel('$|found flux - \mu * flux|$')
    plt.ylabel('N')
    plt.legend(loc='best')
    plt.savefig('/home/ckrawiec/DES/magnification/balrog_N<diffthreshold_flux_mag_'+str(mu)[2:])
    plt.close()

    newsize = new['HALFLIGHTRADIUS_0'] * np.sqrt(mu)
    foundsize = new['HALFLIGHTRADIUS_0_MAG'+str(mu)]
    diffsize = foundsize - newsize
    N = [len(np.where(np.abs(diffsize) < thresh)[0]) for thresh in thresholds]
    plt.plot(thresholds, N)
    plt.xlabel('$found hlr - \mu^{1/2} * hlr$')
    plt.ylabel('N')
    plt.savefig('/home/ckrawiec/DES/magnification/balrog_N<diffthreshold_hlr_mag_'+str(mu)[2:])
    plt.close()

makeplots()
