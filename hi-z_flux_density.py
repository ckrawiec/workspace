#!/usr/bin/env python
"""
Calculate the density of high redshift objects
in flux space of SVA1 GOLD galaxies
using COSMOS photo-z's.
"""
from multiprocessing import Pool
import itertools
import time
import numpy as np
import os
from astropy.io import ascii,fits
import matplotlib.pyplot as plt

num_threads = 4

home_dir = '/home/ckrawiec'
this_file = '{}/git/workspace/hi-z_flux_density.py'.format(home_dir)

#will be overwritten
output_file = '{}/DES/magnification/lbgselect/hi-z_flux_density.fits'.format(home_dir)

sva1_gold_file = '{}/DES/data/sva1_gold_detmodel_gals.fits'.format(home_dir)
sva1_cosmos_file = '{}/DES/data/sva1_coadd_cosmos.fits'.format(home_dir)
cosmos_file = '{}/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'.format(home_dir)

#indices for sva1_gold_detmodel_gals.fits
gold_cosmos15_indices_file = '{}/DES/magnification/lbgselect/gold_cosmos15_indices.txt'.format(home_dir)
gold_no_cosmos_indices_file = '{}/DES/magnification/lbgselect/gold_no_cosmos_indices.txt'.format(home_dir)
cosmos15_indices_file = '{}/DES/magnification/lbgselect/cosmos15_indices.txt'.format(home_dir)


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

def main():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
    print "num_threads="+str(num_threads)
    
    os.environ['OMP_NUM_THREADS']=str(num_threads)

    setup_start = time.time()
    
    #data sets
    sva1_gold = fits.open(sva1_gold_file)[1].data
    sva1_cosmos = fits.open(sva1_cosmos_file)[1].data
    cosmos15 = fits.open(cosmos_file)[1].data

    data_time = time.time()
    print "Loaded data sets in {} s".format(data_time-setup_start)

    gold_cosmos15 = np.loadtxt(gold_cosmos15_indices_file, dtype=int)
    gold_no_cosmos = np.loadtxt(gold_no_cosmos_indices_file, dtype=int)
    cosmos15_indices = np.loadtxt(cosmos15_indices_file, dtype=int)

    #COSMOS2015 photo-z
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z_cosmos = cosmos15['photoz'][cosmos15_indices]
    z0mask = (z_cosmos > 0) & (z_cosmos < 9.9)
    z3mask = (z_cosmos >= 3.) & (z_cosmos < 9.9)
    z4mask = (z_cosmos >= 4.) & (z_cosmos < 9.9)

    #cosmos fluxes and errors from sva1 gold
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

    cosmos_tab_time = time.time()
    hi_z_flux_cosmos = maketable('flux', mask=z4mask, cosmos=True)

    cosmos_tab_end = time.time()
    print "Made COSMOS tables in {}s".format(cosmos_tab_end-cosmos_tab_time)

    hi_z_fluxs = np.array( zip(hi_z_flux_cosmos['g'],
                               hi_z_flux_cosmos['r'],
                               hi_z_flux_cosmos['i'],
                               hi_z_flux_cosmos['z'],
                               hi_z_flux_cosmos['Y']) )

    sva1_tab_time = time.time()
    flux_sva1_gold = maketable('flux')

    sva1_tab_end = time.time()
    print "Made SVA1 table in {}s".format(sva1_tab_end-sva1_tab_time)

    fluxs = np.array( zip(flux_sva1_gold['g'],
                          flux_sva1_gold['r'],
                          flux_sva1_gold['i'],
                          flux_sva1_gold['z'],
                          flux_sva1_gold['Y']) )
    fluxerrs = np.array( zip(flux_sva1_gold['gerr'],
                             flux_sva1_gold['rerr'],
                             flux_sva1_gold['ierr'],
                             flux_sva1_gold['zerr'],
                             flux_sva1_gold['Yerr']) )

    setup_time = time.time()-setup_start
    print "Total setup time took {} s".format(setup_time)
    print " "
    print "# sva1 gold, good region mask, not in COSMOS field: {}".format(len(gold_no_cosmos))
    print "# sva1 gold, good region mask, in COSMOS field: {}".format(len(gold_cosmos15))
    print "# sva1 gold/COSMOS2015 matched, z>4: {}".format(len(gold_cosmos15[z4mask]))
    print "# sva1 gold/COSMOS2015 matched, 0<z<4: {}".format(len(gold_cosmos15[z0mask & ~z4mask]))

    start = time.time()

    N_try = len(fluxs)
    print "Working on {} galaxies...".format(N_try)

    #multiprocessing
    pool = Pool(processes=num_threads)

    n_per_process = int( np.ceil(N_try/num_threads) )
    flux_chunks = [fluxs[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    fluxerr_chunks = [fluxerrs[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]

    hi_results = pool.map(nwrapper, itertools.izip(flux_chunks, fluxerr_chunks, itertools.repeat(hi_z_fluxs)))
    
    hi_final_results = np.concatenate(hi_results)/len(hi_z_fluxs)

    work_time = time.time() - start
    print "Work completed in {} s".format(work_time)

    #write results to fits file
    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(
        [fits.Column(name='coadd_objects_id', format='K', array=sva1_gold['coadd_objects_id'][gold_no_cosmos]),
         fits.Column(name='hi-z_density', format='D', array=hi_final_results)]), nrows=len(hi_final_results))
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Output from {}".format(this_file)
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(output_file, clobber=True)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main()
