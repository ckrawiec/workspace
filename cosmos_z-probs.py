#!/usr/bin/env python  
"""
Calculate the probability that
SVA1 GOLD galaxies belong to redshift
groups using COSMOS photo-z's.
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
this_file = '{}/git/workspace/cosmos_z-probs.py'.format(home_dir)

group = 'lo' #lo (0<z<4) or hi (z>4)
ptype = 'tree' #full or tree
data_type = 'mag' #mag or flux

if ptype='tree':
    pfunc=ptree
elif ptype='full':
    pfunc=p
else:
    raise ValueError('Choose ptype=\'full\' or \'tree\'')

k_near = 20000 #nearest neighbors if ptype=tree

#will be overwritten
output_file = '{}/DES/magnification/lbgselect/{}-z_{}_probability_test.fits'.format(home_dir, group, data_type)

sva1_gold_file = '{}/DES/data/sva1_gold_detmodel_gals.fits'.format(home_dir)
sva1_cosmos_file = '{}/DES/data/sva1_coadd_cosmos.fits'.format(home_dir)
cosmos_file = '{}/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'.format(home_dir)

#indices for sva1_gold_detmodel_gals.fits
gold_cosmos15_indices_file = '{}/DES/magnification/lbgselect/gold_cosmos15_indices.txt'.format(home_dir)
gold_no_cosmos_indices_file = '{}/DES/magnification/lbgselect/gold_no_cosmos_indices.txt'.format(home_dir)
cosmos15_indices_file = '{}/DES/magnification/lbgselect/cosmos15_indices.txt'.format(home_dir)


def pwrapper(args):
    return pfunc(*args)

def p(vals, errs, truevals):
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

def ptree(vals, errs, truevals, knear=k_near):
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
            
            out.append(np.sum(C) * len(truevals)/len(truearr))
            
        return np.array(out)

def main():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
    print "num_threads="+str(num_threads)
    print "Calculating {}-z {} probabilities, using {} integrals".format(group, data_type, ptype)
    
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

    Ntot = len(z_cosmos[z0mask & ~z4mask]) + len(z_cosmos[z4mask])
    
    #cosmos fluxes and errors from sva1 gold
    def maketable(mask=None, cosmos=False, filters=['g','r','i','z','Y']):
        table = {}
        for f in filters:
            if cosmos:
                table[f] = sva1_gold[data_type+'_detmodel_'+f][gold_cosmos15][mask]
                table[f+'err'] = sva1_gold[data_type+'err_detmodel_'+f][gold_cosmos15][mask]
            else:
                table[f] = sva1_gold[datatype+'_detmodel_'+f][gold_no_cosmos]
                table[f+'err'] = sva1_gold[data_type+'err_detmodel_'+f][gold_no_cosmos]
        return table

    if group='hi':
        group_mask=z4mask
    elif group='lo':
        group_mask=(z0mask & ~z4mask)
    else:
        raise ValueError('Unknown group name')

    cosmos_tab_time = time.time()

    data_cosmos = maketable(mask=group_mask, cosmos=True)

    cosmos_tab_end = time.time()
    print "Made COSMOS tables in {}s".format(cosmos_tab_end-cosmos_tab_time)

    truth = np.array( zip(data_cosmos['g'],
                          data_cosmos['r'],
                          data_cosmos['i'],
                          data_cosmos['z'],
                          data_cosmos['Y']) )

    sva1_tab_time = time.time()
    data_sva1_gold = maketable()

    sva1_tab_end = time.time()
    print "Made SVA1 table in {}s".format(sva1_tab_end-sva1_tab_time)

    data = np.array( zip(data_sva1_gold['g'],
                         data_sva1_gold['r'],
                         data_sva1_gold['i'],
                         data_sva1_gold['z'],
                         data_sva1_gold['Y']) )
    data_errs = np.array( zip(data_sva1_gold['gerr'],
                              data_sva1_gold['rerr'],
                              data_sva1_gold['ierr'],
                              data_sva1_gold['zerr'],
                              data_sva1_gold['Yerr']) )

    setup_time = time.time()-setup_start
    print "Total setup time took {} s".format(setup_time)
    print " "
    print "# sva1 gold, good region mask, not in COSMOS field: {}".format(len(gold_no_cosmos))
    print "# sva1 gold, good region mask, in COSMOS field: {}".format(len(gold_cosmos15))
    print "# sva1 gold/COSMOS2015 matched, z>4: {}".format(len(gold_cosmos15[z4mask]))
    print "# sva1 gold/COSMOS2015 matched, 0<z<4: {}".format(len(gold_cosmos15[z0mask & ~z4mask]))
    print "# sva1 gold/COSMOS2015 matched, z>0: {}".format(Ntot)

    start = time.time()

    #multiprocessing
    pool = Pool(processes=num_threads)

    N_try = len(data)
    
    n_per_process = int( np.ceil(N_try/num_threads) )
    data_chunks = [data[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    data_err_chunks = [data_errs[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    
    results = pool.map(pwrapper, itertools.izip(data_chunks, data_err_chunks, itertools.repeat(truth)))
    
    final_results = np.concatenate(results)/Ntot

    work_time = time.time() - start
    print "Work completed in {} s".format(work_time)

    #write results to fits file
    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(
        [fits.Column(name='coadd_objects_id', format='K', array=sva1_gold['coadd_objects_id'][gold_no_cosmos]),
         fits.Column(name=group+'-z_prob', format='D', array=final_results)]), nrows=len(final_results))
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Output from {}".format(this_file)
    if ptype='tree':
        prihdr['NTREE'] = str(k_near)
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(output_file, clobber=True)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main()
