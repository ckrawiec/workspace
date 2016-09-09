#!/usr/bin/env python  
"""
Calculate the probability that
fluxes of input galaxies belong to 
certain redshifts using photo-z's.

usage: ./zprobability.py <data table> <truth table> <output>

<data table> must include the following:
    -SExtractor flux and error outputs

<output> will be a binary fits table
"""
from multiprocessing import Pool
import itertools
import time
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii,fits
from scipy.spatial import ckdtree

num_threads = 4

home_dir = '/home/ckrawiec'
this_file = '{}/git/workspace/cosmos_z-probs.py'.format(home_dir)

ptype = 'full' #full or tree integration
data_type = 'flux' #mag or flux

filters = ['g','r','i','z','Y']

#groups
z_groups = [[0.2,0.8],
            [1.0,3.8],
            [4.0,9.9]]

k_near = 10000 #nearest neighbors if ptype=tree

def pwrapper(args):
    if ptype=='tree':
        return ptree(*args)
    elif ptype=='full':
        return p(*args)
    else:
        raise ValueError('Choose ptype=\'full\' or \'tree\'')
                
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

def main(args):
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
    print "num_threads="+str(num_threads)
    print "Calculating {}-z {} probabilities, using {} integrals".format(group, data_type, ptype)
    
    os.environ['OMP_NUM_THREADS']=str(num_threads)

    setup_start = time.time()
    
    #data
    data = fits.open(args[1])[1].data
    truth_data = fits.open(args[2])[1].data

    data_time = time.time()
    print "Loaded data in {} s".format(data_time-setup_start)

    #COSMOS2015 photo-z
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z = truth['photoz']

    #cosmos fluxes and errors from sva1 gold
    def maketable(mask=None, cosmos=False):
        table = {}
        for f in filters:
            if cosmos:
                table[f] = sva1_gold[data_type+'_detmodel_'+f][gold_cosmos15][mask]
                table[f+'err'] = sva1_gold[data_type+'err_detmodel_'+f][gold_cosmos15][mask]
            else:
                table[f] = sva1_gold[data_type+'_detmodel_'+f][gold_no_cosmos]
                table[f+'err'] = sva1_gold[data_type+'err_detmodel_'+f][gold_no_cosmos]
        return table

    cosmos_tab_time = time.time()

#make table from data
#make table from truth

#for each zpair
#    mask the truth group
#    run the code, save the Ps

#when writing table, make sure to include group definitions
#and info on column names used, input tables, nearest neighbors


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
    N_try = len(data)
    print "Working on {} galaxies...".format(N_try)
    
    pool = Pool(processes=num_threads)

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
    if ptype=='tree':
        prihdr['NTREE'] = str(k_near)
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(output_file)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main()
