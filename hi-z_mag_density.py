"""
Calculate the density of high redshift objects
in magnitude space of SVA1 GOLD galaxies
using COSMOS photo-z's.
"""
import itertools
import time
import esutil
import numpy as np
import os
from astropy.io import ascii,fits
from myutils import match

num_threads = 8

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

    chunks = itertools.izip([vals[i:i+n] for i in xrange(0, len(vals), nchunks)], 
                            [errs[i:i+n] for i in xrange(0, len(vals), nchunks)])
    
    for chunk, errchunk in chunks:
        trueout = np.zeros(len(chunk))
    
        covIs = 1./errchunk**2.
        A = 1./np.sqrt( (2.*np.pi)**len(vals[0])  * np.prod(errchunk**2., axis=1 ) )

        truechunks = (truevals[i:i+ntruechunks] for i in xrange(0, len(truevals), ntruechunks))
        for truechunk in truechunks:
            diff = np.zeros((len(chunk), len(truechunk), len(chunk[0])))
            np.subtract(chunk[:,np.newaxis,:]-truechunk[np.newaxis,:,:], diff)

            B = -0.5 * np.sum(diff**2.*covIs[:,np.newaxis], axis=2)
            C = A[:,np.newaxis] * np.exp(B)
            
            trueout += np.sum(C, axis=1)

        out = np.concatenate((out, trueout))
    
    return out

def maketable(datatype, mask=None, cosmos=False, filters=['g','r','i','z','Y']):
    table = {}
    for f in filters:
        if cosmos:
            table[f] = sva1_gold[datatype+'_detmodel_'+f][gold_cosmos15][mask]
            table[f+'err'] = sva1_gold[datatype+'err_detmodel_'+f][gold_cosmos15][mask]
        else:
            table[f] = sva1_gold[datatype+'_detmodel_'+f][gold_no_cosmos]
            table[f+'err'] = sva1_gold[datatype+'err_detmodel_'+f][gold_no_cosmos]
    return table


def main():
    os.environ['OMP_NUM_THREADS']=str(num_threads)

    setup_start = time.time()
    
    #data sets
    sva1_gold = fits.open('/home/ckrawiec/DES/data/sva1_gold_detmodel_gals.fits')[1].data
    sva1_cosmos = fits.open('/home/ckrawiec/DES/data/sva1_coadd_cosmos.fits')
    cosmos15 = fits.open('/home/ckrawiec/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits')[1].data

    data_time = time.time()
    print "Loaded data sets in {} s".format(data_time-setup_start)

    #good regions mask
    #healpix map
    region_file = '/home/ckrawiec/DES/data/sva1_gold_r1.0_goodregions_04_n4096.fits.gz'
    hpmap = hp.read_map(region_file, nest=True)
    nside = hp.npix2nside(hpmap.size)

    theta = (90.0 - sva1_gold['dec'])*np.pi/180.
    phi = sva1_gold['ra']*np.pi/180.
    pix = hp.ang2pix(nside,theta,phi,nest=True)
    good, = np.where(hpmap[pix] == 1)

    print "Loaded region masks in {} s".format(time.time()-data_time)

    #get gold data for objects in cosmos field
    gold_sort = np.argsort(sva1_gold['coadd_objects_id'][good])
    cosmos_sort = np.argsort(sva1_cosmos['coadd_objects_id'])
    gold_index, cosmos_index = match(sva1_gold['coadd_objects_id'][good][gold_sort],
                                     sva1_cosmos['coadd_objects_id'][cosmos_sort])
    assert all(sva1_gold['coadd_objects_id'][good][gold_sort][gold_index] == sva1_cosmos['coadd_objects_id'][cosmos_sort][cosmos_index])
    gold_cosmos = good[gold_sort][gold_index]
    gold_no_cosmos = good[gold_sort][list(set(np.arange(len(gold_sort))) - set(gold_cosmos))]

    match_time = time.time()

    #match COSMOS field gold objects with COSMOS2015 objects within 1"
    h = esutil.htm.HTM(10)
    h.match(sva1_gold['ra'][gold_cosmos],sva1_gold['dec'][gold_cosmos],
            cosmos15['alpha_j2000'],cosmos15['delta_j2000'],
            radius=1./3600.,
            file='/home/ckrawiec/DES/magnification/lbgselect/match_sva1_gold_cosmos_gals_1arcsec')
    m1 = h.read('/home/ckrawiec/DES/magnification/lbgselect/match_sva1_gold_cosmos_gals_1arcsec')
    gold_m1, cosmos15_m1, merr = zip(*m1)
    gold_cosmos15 = gold_cosmos[gold_m1]

    print "Matched COSMOS objects by ra/dec in {} s".format(time.time()-match_time)

    #COSMOS2015 photo-z
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z_cosmos = cosmos15['photoz'][cosmos15_m1]
    z0mask = (z_cosmos > 0) & (z_cosmos < 9.9)
    z3mask = (z_cosmos >= 3.) & (z_cosmos < 9.9)
    z4mask = (z_cosmos >= 4.) & (z_cosmos < 9.9)

    #cosmos fluxes and errors from sva1 gold

    lo_z_flux_cosmos = maketable('flux', mask=(z0mask & ~z4mask), cosmos=True)
    lo_z_mag_cosmos = maketable('mag', mask=(z0mask & ~z4mask), cosmos=True)

    hi_z_flux_cosmos = maketable('flux', mask=z4mask, cosmos=True)
    hi_z_mag_cosmos = maketable('mag', mask=z4mask, cosmos=True)

    flux_sva1_gold = maketable('flux')
    mag_sva1_gold = maketable('mag')

    mags = np.array( zip(mag_sva1_gold['g'],
                         mag_sva1_gold['r'],
                         mag_sva1_gold['i'],
                         mag_sva1_gold['z']) )
    magerrs = np.array( zip(mag_sva1_gold['gerr'],
                            mag_sva1_gold['rerr'],
                            mag_sva1_gold['ierr'],
                            mag_sva1_gold['zerr']) )
#    fluxes = np.array( zip(flux_sva1_gold['g'],
#                           flux_sva1_gold['r'],
#                           flux_sva1_gold['i'],
#                           flux_sva1_gold['z']) )
#    fluxerrs = np.array( zip(flux_sva1_gold['gerr'],
#                             flux_sva1_gold['rerr'],
#                             flux_sva1_gold['ierr'],
#                             flux_sva1_gold['zerr']) )

#    lo_z_mags = np.array( zip(lo_z_mag_cosmos['g'],
#                            lo_z_mag_cosmos['r'],
#                            lo_z_mag_cosmos['i'],
#                            lo_z_mag_cosmos['z']) )
    hi_z_mags = np.array( zip(hi_z_mag_cosmos['g'],
                            hi_z_mag_cosmos['r'],
                            hi_z_mag_cosmos['i'],
                            hi_z_mag_cosmos['z']) )
#    lo_fluxes = np.array( zip(lo_z_flux_cosmos['g'],
#                            lo_z_flux_cosmos['r'],
#                            lo_z_flux_cosmos['i'],
#                            lo_z_flux_cosmos['z']) )
#    hi_fluxes = np.array( zip(hi_z_flux_cosmos['g'],
#                            hi_z_flux_cosmos['r'],
#                            hi_z_flux_cosmos['i'],
#                            hi_z_flux_cosmos['z']) )

    setup_time = time.time()-setup_start
    print "Total setup time took {} s".format(setup_time)

    start = time.time()
    
    pool = Pool(processes=num_threads)
    
    n_per_process = int( len(gold_no_cosmos)/(num_threads-1) )
    mag_chunks = [mags[i:i+n_per_process] for i in xrange(0, len(mags), n_per_process)])
    magerr_chunks = [magerrs[i:i+n_per_process] for i in xrange(0, len(magerrs), n_per_process)]
    
    results = pool.map(nwrapper, itertools.izip(mag_chunks, magerr_chunks, itertools.repeat(hi_mags)))
    
    final_results = np.concatenate(*results)

    #write results to fits file
    tbhdu = fits.new_table(fits.ColDefs(
        [fits.Column(name='coadd_objects_id', format='K', array=sva1_gold['coadd_objects_id'][gold_no_cosmos]),
         fits.Column(name='hi-z_density', format='D', array=final_results)]))
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Output from /home/ckrawiec/DES/magnification/lbgselect/hi-z_mag_density.py"
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    hdulist.writeto('hi-z_mag_density.fits')
    
    work_time = time.time() - start
    print "Work completed in {} s".format(work_time)

    
if __name__=="__main__":
    main()
