import itertools
import time
import esutil
import numpy as np
from astropy.io import ascii,fits
from myutils import match


def n1(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """
    n = []
    
    for vali, erri in zip(vals, errs):
        covi = np.matrix(np.diag(erri**2.))
        C = 1./np.sqrt( (2.*np.pi)**len(vali)  * np.linalg.det(covi) )
        for item in truevals:
            s = 0
            truearr = np.array(item)
            B = np.exp( -0.5 * np.sum( (np.prod( zip((vali-truearr),
                                                     np.diag(covi.I),
                                                     (vali-truearr))) )) )
            s += C*B
        n.append(s)
        
    return n

def n2(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """

    covIs = 1./errs**2.
    A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errs**2., axis=1 ) )

    B = -0.5 * np.array([np.sum(np.prod(np.array(zip(vals-trueval,covIs,vals-trueval)), axis=1), axis=1)
                         for trueval in truevals])

    C = A*np.exp(B)

    return np.sum(C, axis=0)

def n3(vals, errs, truevals):
    """
    vals, errs, and truevals are lists or arrays of data/error vectors 
    """

    covdets = np.array([np.prod(err**2.) for err in errs])
    covIs = np.array([1./err**2. for err in errs])
    C = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * covdets )
    delvals = np.array([vals-trueval for trueval in truevals])
    B = -0.5 * np.array([np.sum(np.prod(np.array(zip(delval, covIs, delval)), axis=1), axis=1) for delval in delvals])

    D = C*np.exp(B)

    return np.sum(D, axis=0)


def main():
    
    #data sets
    sva1_gold = fits.open('/Users/Christina/COSMOS/sva1_gold_detmodel_gals.fits')[1].data

    #COSMOS2015 photo-z
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z_cosmos = ascii.read('sva1_gold_cosmos_photoz.tab')['redshift']
    z0mask = (z_cosmos > 0) & (z_cosmos < 9.99)
    z3mask = (z_cosmos >= 3.) & (z_cosmos < 9.99)
    z4mask = (z_cosmos >= 4.) & (z_cosmos < 9.99)

    #fluxes and errors from sva1_gold
    flux_cosmos = ascii.read('sva1_gold_cosmos_detmodel_fluxes.tab')
       
    lo_z_flux_cosmos = {'g' : flux_cosmos['g'][z0mask & ~z4mask],
                        'r' : flux_cosmos['r'][z0mask & ~z4mask],
                        'i' : flux_cosmos['i'][z0mask & ~z4mask],
                        'z' : flux_cosmos['z'][z0mask & ~z4mask],
                        'Y' : flux_cosmos['Y'][z0mask & ~z4mask],
                        'gerr' : flux_cosmos['gerr'][z0mask & ~z4mask],
                        'rerr' : flux_cosmos['rerr'][z0mask & ~z4mask],
                        'ierr' : flux_cosmos['ierr'][z0mask & ~z4mask],
                        'zerr' : flux_cosmos['zerr'][z0mask & ~z4mask],
                        'Yerr' : flux_cosmos['Yerr'][z0mask & ~z4mask]}

    hi_z_flux_cosmos = {'g' : flux_cosmos['g'][z4mask],
                        'r' : flux_cosmos['r'][z4mask],
                        'i' : flux_cosmos['i'][z4mask],
                        'z' : flux_cosmos['z'][z4mask],
                        'Y' : flux_cosmos['Y'][z4mask],
                        'gerr' : flux_cosmos['gerr'][z4mask],
                        'rerr' : flux_cosmos['rerr'][z4mask],
                        'ierr' : flux_cosmos['ierr'][z4mask],
                        'zerr' : flux_cosmos['zerr'][z4mask],
                        'Yerr' : flux_cosmos['Yerr'][z4mask]}

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

    flux_sva1_gold = {'g' : sva1_gold['flux_detmodel_g'][good],
                      'r' : sva1_gold['flux_detmodel_r'][good],
                      'i' : sva1_gold['flux_detmodel_i'][good],
                      'z' : sva1_gold['flux_detmodel_z'][good],
                      'Y' : sva1_gold['flux_detmodel_Y'][good],
                      'gerr' : sva1_gold['fluxerr_detmodel_g'][good],
                      'rerr' : sva1_gold['fluxerr_detmodel_r'][good],
                      'ierr' : sva1_gold['fluxerr_detmodel_i'][good],
                      'zerr' : sva1_gold['fluxerr_detmodel_z'][good],
                      'Yerr' : sva1_gold['fluxerr_detmodel_Y'][good]}
    
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

    N_try_lo = 1
    N_try_lo_cosmos = len(lo_z_mag_cosmos['g'])
    N_try_hi = 1
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

    start1lo = time.time()
    n1_lo_z = n1(lomags, lomagerrs, cosmoslomags)
    print "For n1: {}-d mag vector, {} gold galaxies, {} COSMOS lo-z galaxies: {}".format(len(lomags[0]),
                                                                                          N_try_lo,
                                                                                          N_try_lo_cosmos,
                                                                                          time.time()-start1lo)
    start2lo = time.time()
    n2_lo_z = n2(lomags, lomagerrs, cosmoslomags)
    print "For n2: {}-d mag vector, {} gold galaxies, {} COSMOS lo-z galaxies: {}".format(len(lomags[0]),
                                                                                          N_try_lo,
                                                                                          N_try_lo_cosmos,
                                                                                          time.time()-start2lo)
    start3lo = time.time()
    n3_lo_z = n3(lomags, lomagerrs, cosmoslomags)
    print "For n3: {}-d mag vector, {} gold galaxies, {} COSMOS lo-z galaxies: {}".format(len(lomags[0]),
                                                                                          N_try_lo,
                                                                                          N_try_lo_cosmos,
                                                                                          time.time()-start3lo)
    start1hi = time.time()
    n1_hi_z = n1(himags, himagerrs, cosmoshimags)
    print "For n1: {}-d mag vector, {} gold galaxies, {} COSMOS hi-z galaxies: {}".format(len(himags[0]),
                                                                                          N_try_hi,
                                                                                          N_try_hi_cosmos,
                                                                                          time.time()-start1hi)
    start2hi = time.time()
    n2_lo_z = n2(himags, himagerrs, cosmoshimags)
    print "For n2: {}-d mag vector, {} gold galaxies, {} COSMOS hi-z galaxies: {}".format(len(himags[0]),
                                                                                          N_try_hi,
                                                                                          N_try_hi_cosmos,
                                                                                          time.time()-start2hi)
    start3hi = time.time()
    n3_lo_z = n3(himags, himagerrs, cosmoshimags)
    print "For n3: {}-d mag vector, {} gold galaxies, {} COSMOS hi-z galaxies: {}".format(len(himags[0]),
                                                                                          N_try_hi,
                                                                                          N_try_hi_cosmos,
                                                                                          time.time()-start2hi)
    
    
if __name__=="__main__":
    main()
