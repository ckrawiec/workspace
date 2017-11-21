import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import ConfigParser
from scipy.spatial import ckdtree
from astropy.io import fits
from astropy.table import Table

def parseconfig(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    params['source_file'] = config.get('I/O','source_file')
    params['lens_file']   = config.get('I/O','lens_file')
    params['random_file'] = config.get('I/O','random_file')

    params['source_ra']   = config.get('Columns', 'source_ra')
    params['source_dec']  = config.get('Columns', 'source_dec')
    params['lens_ra']     = config.get('Columns', 'lens_ra')
    params['lens_dec']    = config.get('Columns', 'lens_dec')
    params['rand_ra']     = config.get('Columns', 'random_ra')
    params['rand_dec']    = config.get('Columns', 'random_dec')

    return params

def main():
    params = parseconfig('findpairs.config')
    
    #radii - same units as positions
    r_range = np.logspace(-2, 0, 5)

    #read files
    start_d = time.time()
    source_dat = fits.open(params['source_file'])[1].data
    lens_dat = fits.open(params['lens_file'])[1].data
    rand_dat = fits.open(params['random_file'])[1].data
    end_d = time.time()

    sys.stderr.write('Data read in {}s\n'.format(end_d-start_d))
    
    sys.stderr.write('Sources: {}\nLenses: {}\n'.format(len(source_dat), len(lens_dat)))
    
    #positions
    source_pos = np.array( zip(source_dat[params['source_ra']],
                               source_dat[params['source_dec']]) )
    lens_pos = np.array( zip(lens_dat[params['lens_ra']],
                             lens_dat[params['lens_dec']]) )
    rand_pos = np.array( zip(rand_dat[params['rand_ra']],
                             rand_dat[params['rand_dec']]) )
    
    #make tree
    start_tree = time.time()
    source_tree = ckdtree.cKDTree(source_pos)
    rand_tree = ckdtree.cKDTree(rand_pos)
    end_tree = time.time()
    
    sys.stderr.write('Trees created in {}s\n'.format(end_tree-start_tree))
    
    start_q = time.time()
    
    sys.stderr.write('Starting queries...\n')
    
    pair_dict = {}
    rpair_dict = {}
    
    #for each radius, query for all sources around lenses
    for ri in range(len(r_range)):
        pair_dict[r_range[ri]] = []
        rpair_dict[r_range[ri]] = []
        sys.stderr.write('    working on r={}\n'.format(r_range[ri]))
        pairs2 = source_tree.query_ball_point(lens_pos, r=r_range[ri])
        rpairs2 = rand_tree.query_ball_point(lens_pos, r=r_range[ri])
        if ri==0:
            pairs = len(np.hstack(pairs2))
            rpairs = len(np.hstack(rpairs2))
        else:
            pairs1 = source_tree.query_ball_point(lens_pos, r=r_range[ri-1])
            rpairs1 = rand_tree.query_ball_point(lens_pos, r=r_range[ri-1])
            pairs = len(np.hstack(pairs2))-len(np.hstack(pairs1))
            rpairs = len(np.hstack(rpairs2))-len(np.hstack(rpairs1))
            
        pair_dict[r_range[ri]].append(pairs)
        rpair_dict[r_range[ri]].append(rpairs)
    
    #randoms & MASKS
    
    end_q = time.time()
    
    sys.stderr.write('Time for queries: {}s\n'.format(end_q-start_q))
    
    for k in pair_dict.keys():
        print 'r={}: {} pairs'.format(k, pair_dict[k])
        print '      {} random pairs'.format(rpair_dict[k])
    
    plt.scatter(pair_dict.keys(), pair_dict.values(), c='b', label='gold')
    plt.scatter(rpair_dict.keys(), rpair_dict.values(), c='g', label='balrog')
    
    plt.xlabel('r (deg)')
    plt.ylabel('source-lens pairs')
    plt.xscale('log')
    #plt.yscale('log')
    plt.show()

    
if __name__=="__main__":
    main()
