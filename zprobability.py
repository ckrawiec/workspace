usage = """
Calculate the Bayesian probability that
fluxes of input galaxies belong to 
certain redshift ranges using photo-z's.

usage: python zprobability.py <data table> <truth table> <output> <start index> <end index>

<data table> must include the following:
    -SExtractor flux and error outputs
<truth table> must include the following:
    -SExtractor flux and error outputs
    -column named 'photoz'

<output> will be a binary fits table and will overwrite existing file

<start/end index> are optional indices to only work on a slice of the data table
"""
import itertools
import time
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from astropy.io import ascii,fits
from scipy.spatial import ckdtree

num_threads = 8

ptype = 'tree' #full or tree integration
data_type = 'FLUX' #MAG or FLUX

filters = ['G','R','I','Z','Y']

#groups
z_groups = [[0.01,4.0],
            [4.0,9.9]]

k_near = 10000 #nearest neighbors if ptype=tree

def pwrapper(args):
    if ptype=='tree':
        #if truth array shorter than k_near, just do full integration
        if len(args[-1]) > k_near:
            return ptree(*args)
        else:
            print "Length of truth array <= {}, using brute force instead of tree".format(k_near)
            return p(*args)
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
    if len(args)<4:
        print "Missing command line argument."
        print usage
        sys.exit()

    data_file = args[1]
    template_file = args[2]
    output_file = args[3]
    if len(args)==6:
        start_index = int(args[4])
        end_index = int(args[5])
    else:
        start_index = 0
        end_index = None

    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
    print "num_threads="+str(num_threads)
    print "Calculating {} probabilities, using {} integrals".format(data_type, ptype)
    
    os.environ['OMP_NUM_THREADS']=str(num_threads)

    setup_start = time.time()
    
    #open tables
    data = fits.open(data_file)[1].data
    templates = fits.open(template_file)[1].data

    data_time = time.time()
    print "Loaded data in {} s".format(data_time-setup_start)

    #redshifts
    z = templates['photoz']

    #data vectors
    data_zip = np.array( zip( *[data[data_type+'_detmodel_'+f][start_index:end_index] for f in filters] ) )
    err_zip = np.array( zip( *[data[data_type+'err_detmodel_'+f][start_index:end_index] for f in filters] ) )

    data_ids = data['COADD_OBJECTS_ID'][start_index:end_index]
    del data

    setup_time = time.time()-setup_start
    print "Total setup time took {} s".format(setup_time)
    
    P_dict = {}

    N_try = len(data_zip)
    print "Working on {} galaxies (table indices {}-{}) ...".format(N_try, start_index, end_index)

    n_per_process = int( np.ceil(N_try/num_threads) )
    data_chunks = [data_zip[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    err_chunks = [err_zip[i:i+n_per_process] for i in xrange(0, N_try, n_per_process)]
    
    del data_zip
    del err_zip
    
    #multiprocessing
    pool = Pool(processes=num_threads)

    for z_group in z_groups:
        z_mask = (z >= np.min(z_group)) & (z < np.max(z_group))

        template_zip = np.array( zip( *[templates[data_type+'_detmodel_'+f][z_mask] for f in filters] ) )
        print "# of COSMOS templates in z group {}: {}".format(str(z_group), len(template_zip))
    
        start = time.time()
        
        results = pool.map(pwrapper, itertools.izip(data_chunks, err_chunks, itertools.repeat(template_zip)))
        
        P_dict[str(z_group)] = np.concatenate(results)

        work_time = time.time() - start
        print "Work completed in {} s".format(work_time)

    pool.close()
        
    #write results to fits file
    col_defs = [fits.Column(name='COADD_OBJECTS_ID', format='K', array=data_ids)]

    P_norm = np.zeros(N_try)
    for k in P_dict.keys():
        P_norm += P_dict[k]
    for z_group in z_groups:
        col_defs.append(fits.Column(name='P'+str(z_group), format='D', array=P_dict[str(z_group)]/P_norm))

    col_defs.append(fits.Column(name='PNORM', format='D', array=P_norm))

    pri_hdr = fits.Header()
    tb_hdr = fits.Header()
    tb_hdr['COMMENT'] = "Bayesian redshift probabilities for data in {} using photo-zs of templates from {}. Data vectors/errors were comprised of {}(ERR)_DETMODEL_% for % in {}. Columns reported here are \'P[zmin, zmax]\'".format(data_file, template_file, data_type, filters)
    if ptype=='tree':
        tb_hdr['NTREE'] = str(k_near)

    pri_hdu = fits.PrimaryHDU(header=pri_hdr)
    tb_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(col_defs), nrows=N_try, header=tb_hdr)
    hdu_list = fits.HDUList([pri_hdu, tb_hdu])
    hdu_list.writeto(output_file, clobber=True)
    
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now
        
if __name__=="__main__":
    main(sys.argv)
