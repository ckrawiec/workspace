import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import ckdtree
from astropy.table import Table

num_threads = 1

def printtime():
    now = time.strftime("%Y-%m-%d %H:%M")
    print "#"+now

printtime()

balrog_file = sys.argv[1]
output = sys.argv[2]
mu = float(sys.argv[3])

balrog = Table.read(balrog_file)

flux = np.array([balrog['FLUX_NOISELESS_G'],
                 balrog['FLUX_NOISELESS_R'],
                 balrog['FLUX_NOISELESS_I'],
                 balrog['FLUX_NOISELESS_Z'],
                 balrog['FLUX_NOISELESS_Y']])

n_try = len(balrog)

mag_flux = flux[:n_try] * mu
mag_size = balrog['HALFLIGHTRADIUS_0'][:n_try] * np.sqrt(mu)

mag_data = np.array( zip(*np.vstack([mag_flux, mag_size])) )
data = np.array( zip(*np.vstack([flux, balrog['HALFLIGHTRADIUS_0']])) )

del flux
del mag_flux

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

del data_chunks

print "query time: {}s".format(map_end-map_start)

mag_ids = np.concatenate(results)

print "first element of magnified data: {}, matched vector id: {}".format(mag_data[0], mag_ids[0])

t = Table()
t['BALROG_INDEX'] = balrog['BALROG_INDEX'][:n_try]

indices = [i[1] for i in mag_ids]
t['INDEX_MAG'+str(mu)] = t['BALROG_INDEX'][indices]
t['FLUX_NOISELESS_G_MAG'+str(mu)] = balrog['FLUX_NOISELESS_G'][indices]
t['FLUX_NOISELESS_R_MAG'+str(mu)] = balrog['FLUX_NOISELESS_R'][indices]
t['FLUX_NOISELESS_I_MAG'+str(mu)] = balrog['FLUX_NOISELESS_I'][indices]
t['FLUX_NOISELESS_Z_MAG'+str(mu)] = balrog['FLUX_NOISELESS_Z'][indices]
t['FLUX_NOISELESS_Y_MAG'+str(mu)] = balrog['FLUX_NOISELESS_Y'][indices]
t['HALFLIGHTRADIUS_0_MAG'+str(mu)] = balrog['HALFLIGHTRADIUS_0'][indices]
t['D_MAG'+str(mu)] = [j[0] for j in mag_ids]

if os.path.exists(output):
    spl  = os.path.splitext(output)
    output = "{}-new{}".format(spl[0],spl[1])
    print "output file exists: output will be written to {}".format(output)
    
t.write(output)
print "output written to {}".format(output)

printtime()
