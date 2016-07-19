import itertools
import time
import esutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ckdtree
from astropy.io import ascii,fits
from myutils import match

def ntree(vals, errs, truevals):
    knear = 1500
    truetree = ckdtree.cKDTree(truevals)

    out = []
    for val, err in zip(vals, errs):
        dnear, inear = truetree.query(val, k=knear)
        truearr = truetree.data[inear]

        covI = 1./err**2.
        A = 1./np.sqrt( (2.*np.pi)**len(val)  * np.prod(err**2.) )

        diff = val-truearr
        B = -0.5 * np.sum(diff**2.*covI, axis=1)
        C = A * np.exp(B)
        
        out.append(np.sum(C))

    return np.array(out)


def nwrapper(args):
    return n(*args)

def n(vals, errs, truevals):
    """
    sum the gaussian likelihoods L(vals|truevals) over truevals using the errs on vals
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """
    out = np.array([])

    nchunks = 500
    ntruechunks = 1000

    chunks = itertools.izip([vals[i:i+nchunks] for i in xrange(0, len(vals), nchunks)], 
                            [errs[i:i+nchunks] for i in xrange(0, len(vals), nchunks)])
    
    for chunk, errchunk in chunks:
        trueout = np.zeros(len(chunk))
    
        covIs = 1./errchunk**2.
        A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errchunk**2., axis=1 ) )

        truechunks = (truevals[i:i+ntruechunks] for i in xrange(0, len(truevals), ntruechunks))
        for truechunk in truechunks:
            diff = chunk[:,np.newaxis,:]-truechunk[np.newaxis,:,:]

            B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)
            C = A[:,np.newaxis] * np.exp(B)
            
            trueout += np.sum(C, axis=1)

        out = np.concatenate((out, trueout))
    
    return out


def n1(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """
    n = 500
    chunks = itertools.izip([vals[i:i+n] for i in xrange(0, len(vals), n)], 
                            [errs[i:i+n] for i in xrange(0, len(vals), n)])
    
    ntrue = 1000

    out = np.array([])
    
    for chunk, errchunk in chunks:
        trueout = np.zeros(len(chunk))
    
        covIs = 1./errchunk**2.
        A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errchunk**2., axis=1 ) )

        truechunks = (truevals[i:i+ntrue] for i in xrange(0, len(truevals), ntrue))
        for truechunk in truechunks:
            diff = chunk[:,np.newaxis,:]-truechunk[np.newaxis,:,:]
            B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)
            
            C = A[:,np.newaxis] * np.exp(B)
            trueout += np.sum(C, axis=1)

        out = np.concatenate((out, trueout))

    return out

def n2(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """
    covIs = 1./errs**2.
    A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errs**2., axis=1 ) )

    out = np.zeros(len(vals))
    
    n = 1000
    truechunks = (truevals[i:i+n] for i in xrange(0, len(truevals), n))
    for truechunk in truechunks:
        diff = np.zeros((len(vals), len(truechunk),  len(vals[0])))
        np.subtract(vals[:,np.newaxis,:], truechunk[np.newaxis,:,:], diff)
        B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)

        C = A[:,np.newaxis] * np.exp(B)
        out += np.sum(C, axis=1)

    return out

def n3(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """
    covIs = 1./errs**2.
    A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errs**2., axis=1 ) )

    out = np.zeros(len(vals))
    
    n = 1000
    truechunks = (truevals[i:i+n] for i in xrange(0, len(truevals), n))
    for truechunk in truechunks:
        diff = vals[:,np.newaxis,:]-truechunk[np.newaxis,:,:]
        B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)

        C = A[:,np.newaxis] * np.exp(B)
        out += np.sum(C, axis=1)

    return out

def main():
    
    #data sets
    sva1_gold = fits.open('/Users/Christina/DES/data/sva1_gold_detmodel_gals.fits')[1].data

    #COSMOS2015 photo-z
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z_cosmos = ascii.read('sva1_gold_cosmos_photoz.tab')['redshift']
    z0mask = (z_cosmos > 0) & (z_cosmos < 9.9)
    z3mask = (z_cosmos >= 3.) & (z_cosmos < 9.9)
    z4mask = (z_cosmos >= 4.) & (z_cosmos < 9.9)

    #fluxes and errors from sva1_gold
    flux_cosmos = ascii.read('sva1_gold_cosmos_detmodel_fluxes.tab')
       
#    lo_z_flux_cosmos = {'g' : flux_cosmos['g'][z0mask & ~z4mask],
#                        'r' : flux_cosmos['r'][z0mask & ~z4mask],
#                        'i' : flux_cosmos['i'][z0mask & ~z4mask],
#                        'z' : flux_cosmos['z'][z0mask & ~z4mask],
#                        'Y' : flux_cosmos['Y'][z0mask & ~z4mask],
#                        'gerr' : flux_cosmos['gerr'][z0mask & ~z4mask],
#                        'rerr' : flux_cosmos['rerr'][z0mask & ~z4mask],
#                        'ierr' : flux_cosmos['ierr'][z0mask & ~z4mask],
#                        'zerr' : flux_cosmos['zerr'][z0mask & ~z4mask],
#                        'Yerr' : flux_cosmos['Yerr'][z0mask & ~z4mask]}

#    hi_z_flux_cosmos = {'g' : flux_cosmos['g'][z4mask],
#                        'r' : flux_cosmos['r'][z4mask],
#                        'i' : flux_cosmos['i'][z4mask],
#                        'z' : flux_cosmos['z'][z4mask],
#                        'Y' : flux_cosmos['Y'][z4mask],
#                        'gerr' : flux_cosmos['gerr'][z4mask],
#                        'rerr' : flux_cosmos['rerr'][z4mask],
#                        'ierr' : flux_cosmos['ierr'][z4mask],
#                        'zerr' : flux_cosmos['zerr'][z4mask],
#                        'Yerr' : flux_cosmos['Yerr'][z4mask]}

    mag_cosmos = ascii.read('sva1_gold_cosmos_detmodel_mags.tab')
                
    lo_z_mag_cosmos = {'g' : mag_cosmos['g'][z0mask & ~z4mask],
                       'r' : mag_cosmos['r'][z0mask & ~z4mask],
                       'i' : mag_cosmos['i'][z0mask & ~z4mask],
                       'z' : mag_cosmos['z'][z0mask & ~z4mask],
                       'Y' : mag_cosmos['Y'][z0mask & ~z4mask],
                       'gerr' : mag_cosmos['gerr'][z0mask & ~z4mask],
                       'rerr' : mag_cosmos['rerr'][z0mask & ~z4mask],
                       'ierr' : mag_cosmos['ierr'][z0mask & ~z4mask],
                       'zerr' : mag_cosmos['zerr'][z0mask & ~z4mask],
                       'Yerr' : mag_cosmos['Yerr'][z0mask & ~z4mask]}

    hi_z_mag_cosmos = {'g' : mag_cosmos['g'][z4mask],
                       'r' : mag_cosmos['r'][z4mask],
                       'i' : mag_cosmos['i'][z4mask],
                       'z' : mag_cosmos['z'][z4mask],
                       'Y' : mag_cosmos['Y'][z4mask],
                       'gerr' : mag_cosmos['gerr'][z4mask],
                       'rerr' : mag_cosmos['rerr'][z4mask],
                       'ierr' : mag_cosmos['ierr'][z4mask],
                       'zerr' : mag_cosmos['zerr'][z4mask],
                       'Yerr' : mag_cosmos['Yerr'][z4mask]}

    good = np.loadtxt('sva1_gold_goodregions_indices.txt.gz', dtype=int)

#    flux_sva1_gold = {'g' : sva1_gold['flux_detmodel_g'][good],
#                      'r' : sva1_gold['flux_detmodel_r'][good],
#                      'i' : sva1_gold['flux_detmodel_i'][good],
#                      'z' : sva1_gold['flux_detmodel_z'][good],
#                      'Y' : sva1_gold['flux_detmodel_Y'][good],
#                      'gerr' : sva1_gold['fluxerr_detmodel_g'][good],
#                      'rerr' : sva1_gold['fluxerr_detmodel_r'][good],
#                      'ierr' : sva1_gold['fluxerr_detmodel_i'][good],
#                      'zerr' : sva1_gold['fluxerr_detmodel_z'][good],
#                      'Yerr' : sva1_gold['fluxerr_detmodel_Y'][good]}
    
    mag_sva1_gold = {'g' : sva1_gold['mag_detmodel_g'][good],
                     'r' : sva1_gold['mag_detmodel_r'][good],
                     'i' : sva1_gold['mag_detmodel_i'][good],
                     'z' : sva1_gold['mag_detmodel_z'][good],
                     'Y' : sva1_gold['mag_detmodel_Y'][good],
                     'gerr' : sva1_gold['magerr_detmodel_g'][good],
                     'rerr' : sva1_gold['magerr_detmodel_r'][good],
                     'ierr' : sva1_gold['magerr_detmodel_i'][good],
                     'zerr' : sva1_gold['magerr_detmodel_z'][good],
                     'Yerr' : sva1_gold['magerr_detmodel_Y'][good]}

    N_try_lo = 1000
    N_try_lo_cosmos = len(lo_z_mag_cosmos['g'])
    N_try_hi = 1000
    N_try_hi_cosmos = len(hi_z_mag_cosmos['g'])

    lomags = np.array( zip(mag_sva1_gold['g'][:N_try_lo],
                           mag_sva1_gold['r'][:N_try_lo],
                           mag_sva1_gold['i'][:N_try_lo]) )
    lomagerrs = np.array( zip(mag_sva1_gold['gerr'][:N_try_lo],
                              mag_sva1_gold['rerr'][:N_try_lo],
                              mag_sva1_gold['ierr'][:N_try_lo]) )
    cosmoslomags = np.array( zip(lo_z_mag_cosmos['g'][:N_try_lo_cosmos],
                                 lo_z_mag_cosmos['r'][:N_try_lo_cosmos],
                                 lo_z_mag_cosmos['i'][:N_try_lo_cosmos]) )

    himags = np.array( zip(mag_sva1_gold['g'][:N_try_hi],
                           mag_sva1_gold['r'][:N_try_hi],
                           mag_sva1_gold['i'][:N_try_hi]) )
    himagerrs = np.array( zip(mag_sva1_gold['gerr'][:N_try_hi],
                              mag_sva1_gold['rerr'][:N_try_hi],
                              mag_sva1_gold['ierr'][:N_try_hi]) )
    cosmoshimags = np.array( zip(hi_z_mag_cosmos['g'][:N_try_hi_cosmos],
                                 hi_z_mag_cosmos['r'][:N_try_hi_cosmos],
                                 hi_z_mag_cosmos['i'][:N_try_hi_cosmos]) )


    def ntest(vals, errs, truevals, nfunc, N_try, N_try_cosmos):
        start = time.time()
        out = nfunc(vals, errs, truevals)
        end = time.time()
        print "For {}: {}-d mag vector, {} gold galaxies, {} COSMOS galaxies: {}".format(nfunc.func_name,
                                                                                         len(vals[0]),
                                                                                         N_try,
                                                                                         N_try_cosmos,
                                                                                         end-start)

        return out

    def plotn(n_lo, n_hi, save_file):
        plt.scatter( mag_sva1_gold['r'][:N_try_hi] - mag_sva1_gold['i'][:N_try_hi],
                     mag_sva1_gold['g'][:N_try_hi] - mag_sva1_gold['r'][:N_try_hi],
                     c=(n_hi/N_try_hi_cosmos)/( (n_hi/N_try_hi_cosmos)+(n_lo/N_try_lo_cosmos) ),
                     edgecolor='none')
        plt.xlabel('mag_detmodel_r - mag_detmodel_i')
        plt.ylabel('mag_detmodel_g - mag_detmodel_r')
        plt.colorbar(label='hi/(hi+lo)')
        plt.xlim(-1,4)
        plt.ylim(-1,4)
        plt.savefig(save_file)
        plt.close()

    n_lo_z = ntest(lomags, lomagerrs, cosmoslomags, n, N_try_lo, N_try_lo_cosmos)
    n1_lo_z = ntest(lomags, lomagerrs, cosmoslomags, n1, N_try_lo, N_try_lo_cosmos)
    n2_lo_z = ntest(lomags, lomagerrs, cosmoslomags, n2, N_try_lo, N_try_lo_cosmos)
    n3_lo_z = ntest(lomags, lomagerrs, cosmoslomags, n3, N_try_lo, N_try_lo_cosmos)
    ntree_lo_z = ntest(lomags, lomagerrs, cosmoslomags, ntree, N_try_lo, N_try_lo_cosmos)

    n_hi_z = ntest(himags, himagerrs, cosmoshimags, n, N_try_hi, N_try_hi_cosmos)
    n1_hi_z = ntest(himags, himagerrs, cosmoshimags, n1, N_try_hi, N_try_hi_cosmos)
    n2_hi_z = ntest(himags, himagerrs, cosmoshimags, n2, N_try_hi, N_try_hi_cosmos)
    n3_hi_z = ntest(himags, himagerrs, cosmoshimags, n3, N_try_hi, N_try_hi_cosmos)
    ntree_hi_z = ntest(himags, himagerrs, cosmoshimags, ntree, N_try_hi, N_try_hi_cosmos)


    plotn(n_lo_z, n_hi_z, 'hi_z_cosmos_mag_densities_ri_gr_n')
    plotn(n1_lo_z, n1_hi_z, 'hi_z_cosmos_mag_densities_ri_gr_n1')
    plotn(n2_lo_z, n2_hi_z, 'hi_z_cosmos_mag_densities_ri_gr_n2')
    plotn(n3_lo_z, n3_hi_z, 'hi_z_cosmos_mag_densities_ri_gr_n3')
    plotn(ntree_lo_z, ntree_hi_z, 'hi_z_cosmos_mag_densities_ri_gr_ntree')

    if (np.all(n1_lo_z==n2_lo_z) and np.all(n2_lo_z==n3_lo_z) and np.all(n1_hi_z==n_hi_z)):
        check = 'True'
    else:
        check = 'False'

    print "Outputs are the same for n,n1,n2,n3?:" + check
    

    
if __name__=="__main__":
    main()
