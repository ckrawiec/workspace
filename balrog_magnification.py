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
mean_size = np.mean(mag_size)
sigma_size = np.std(mag_size)

new_mag_flux = np.array([(mfi-np.mean(mfi))/np.std(mfi) for mfi in mag_flux])
new_flux = np.array([(fi-np.mean(mfi))/np.std(mfi) for fi,mfi in zip(flux,mag_flux)])
new_mag_size = (mag_size - mean_size) / sigma_size
new_size = (balrog['HALFLIGHTRADIUS_0'] - mean_size) / sigma_size

del mag_flux
del mag_size

mag_data = np.array( zip(*np.vstack([new_mag_flux, new_mag_size])) )
data = np.array( zip(*np.vstack([new_flux, new_size])) )

del flux

g = balrog['FLUX_NOISELESS_G']
r = balrog['FLUX_NOISELESS_R']
i = balrog['FLUX_NOISELESS_I']
z = balrog['FLUX_NOISELESS_Z']
y = balrog['FLUX_NOISELESS_Y']
size = balrog['HALFLIGHTRADIUS_0']
index = balrog['BALROG_INDEX']

del balrog

st = time.time()
tree = ckdtree.cKDTree(data[:n_try], balanced_tree=False)
en = time.time()

del data

print "tree build time ({} elements): {}s".format(n_try, en-st)

def savequery(data_chunk):
    return tree.query(data_chunk)

#multiprocessing
n_per_process = int( np.ceil(n_try/num_threads) )

data_chunks = (mag_data[i:i+n_per_process] for i in xrange(0, n_try, n_per_process))

print "before pool"
pool = Pool(processes=num_threads)
print "after pool"

map_start = time.time()
results = pool.map(savequery, data_chunks)
map_end = time.time()

pool.close()

del data_chunks

print "query time: {}s".format(map_end-map_start)

dists, indices = np.hstack(results)

print "first element of magnified data: {}, matched vector id: {}".format(mag_data[0], indices[0])

del mag_data

t = Table()
t['BALROG_INDEX'] = index[:n_try]

indices = [int(ii) for ii in indices]

t['BALROG_INDEX_MAG'+str(mu)] = t['BALROG_INDEX'][indices]
t['FLUX_NOISELESS_G_MAG'+str(mu)] = g[indices]
t['FLUX_NOISELESS_R_MAG'+str(mu)] = r[indices]
t['FLUX_NOISELESS_I_MAG'+str(mu)] = i[indices]
t['FLUX_NOISELESS_Z_MAG'+str(mu)] = z[indices]
t['FLUX_NOISELESS_Y_MAG'+str(mu)] = y[indices]
t['HALFLIGHTRADIUS_0_MAG'+str(mu)] = size[indices]
t['D_MAG'+str(mu)] = dists
t['I_MAG'+str(mu)] = indices

if os.path.exists(output):
    spl  = os.path.splitext(output)
    output = "{}-new{}".format(spl[0],spl[1])
    print "output file exists: output will be written to {}".format(output)
    
t.write(output)
print "output written to {}".format(output)

printtime()
