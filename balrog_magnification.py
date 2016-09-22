import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import ckdtree
from astropy.table import Table

num_threads = 2

def printtime():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now

printtime()

balrog_file = sys.argv[1]
output = sys.argv[2]
mu = sys.argv[3]

balrog = Table.read(balrog_file)

data = np.array( zip(balrog['FLUX_NOISELESS_G'],
                     balrog['FLUX_NOISELESS_R'],
                     balrog['FLUX_NOISELESS_I'],
                     balrog['FLUX_NOISELESS_Z'],
                     balrog['FLUX_NOISELESS_Y']) )

n_try = len(data)

mag_data = data[:n_try] * mu

st = time.time()
tree = ckdtree.cKDTree(data[:n_try], balanced_tree=False)
en = time.time()

print "tree build time ({} elements): {}s".format(n_try, en-st)

def savequery(data_chunk):
    ids = []
    for d in data_chunk:
        ids.append(tree.query(d))
    return ids

#multiprocessing
n_per_process = int( np.ceil(n_try/num_threads) )

data_chunks = (mag_data[i:i+n_per_process] for i in xrange(0, n_try, n_per_process))

print "before pool"
pool = Pool(processes=num_threads)
print "after pool"

map_start = time.time()
results = pool.map(savequery, data_chunks)
map_end = time.time()

print "query time: {}s".format(map_end-map_start)

mag_ids = np.concatenate(results)

print "first element of magnified data: {}, matched vector id: {}".format(mag_data[0], mag_ids[0])

t = Table()
t['BALROG_INDEX'] = balrog['BALROG_INDEX'][:n_try]
t['ID_MAG'+str(mu)] = mag_ids

if os.path.exists(output):
    spl  = os.path.splitext(output)
    output = "{}-new{}".format(spl[0],spl[1])
    print "output file exists: output will be written to {}".format(output)
    
t.write(output)
print "output written to {}".format(output)

printtime()
