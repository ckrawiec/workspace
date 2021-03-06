#!/usr/bin/env python
"""
Calculate the density of low redshift objects
in flux space of SVA1 GOLD galaxies
using COSMOS photo-z's.
"""
import itertools
import time
import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import ckdtree
from astropy.io import ascii,fits

num_threads = 4

#how many nearest neighbors to use in likelihood sum
k_near = 10000

home_dir = '/home/ckrawiec'
this_file = '{}/git/workspace/lo-z_flux_density_tree.py'.format(home_dir)

#will be overwritten
output_file = '{}/DES/magnification/lbgselect/lo-z_flux_density_tree.fits'.format(home_dir)

sva1_gold_file = '{}/DES/data/sva1_gold_detmodel_gals.fits'.format(home_dir)
sva1_cosmos_file = '{}/DES/data/sva1_coadd_cosmos.fits'.format(home_dir)
cosmos_file = '{}/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'.format(home_dir)

#indices for sva1_gold_detmodel_gals.fits
gold_cosmos15_indices_file = '{}/DES/magnification/lbgselect/gold_cosmos15_indices.txt'.format(home_dir)
gold_no_cosmos_indices_file = '{}/DES/magnification/lbgselect/gold_no_cosmos_indices.txt'.format(home_dir)
cosmos15_indices_file = '{}/DES/magnification/lbgselect/cosmos15_indices.txt'.format(home_dir)


def nwrapper(args):
    return ntree(*args)

def ntree(vals, errs, truevals, knear=k_near):
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
    lo_z_flux_cosmos = maketable('flux', mask=(z0mask & ~z4mask), cosmos=True)

    cosmos_tab_end = time.time()
    print "Made COSMOS tables in {}s".format(cosmos_tab_end-cosmos_tab_time)

    lo_z_fluxs = np.array( zip(lo_z_flux_cosmos['g'],
                               lo_z_flux_cosmos['r'],
                               lo_z_flux_cosmos['i'],
                               lo_z_flux_cosmos['z'],
                               lo_z_flux_cosmos['Y']) )

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

    #multiprocessing
    pool = Pool(processes=num_threads)

    N_try = len(fluxs)

    n_per_process = int( np.ceil(N_try/num_threads) )
    flux_chunks = [fluxs[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    fluxerr_chunks = [fluxerrs[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]

    lo_results = pool.map(nwrapper, itertools.izip(flux_chunks, fluxerr_chunks, itertools.repeat(lo_z_fluxs)))
    
    lo_final_results = np.concatenate(lo_results)/len(lo_z_fluxs)

    work_time = time.time() - start
    print "Work completed in {} s".format(work_time)

    #write results to fits file
    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(
        [fits.Column(name='coadd_objects_id', format='K', array=sva1_gold['coadd_objects_id'][gold_no_cosmos]),
         fits.Column(name='lo-z_density', format='D', array=lo_final_results)]), nrows=len(lo_final_results))
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Output from {}".format(this_file)
    prihdr['NNEIGHBORS'] = (str(k_near), 'number of nearest neighbors used by kd-tree')
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(output_file, clobber=True)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main()
