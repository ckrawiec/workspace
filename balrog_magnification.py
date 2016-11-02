import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ckdtree
from astropy.table import Table

def printtime():
    now = time.strftime("%Y-%m-%d %H:%M")
    sys.stderr.write('#'+now+'\n')

printtime()

balrog_file = sys.argv[1]
output = sys.argv[2]
mu = float(sys.argv[3])

balrog = Table.read(balrog_file)

sys.stderr.write('table read\n')

n_try = 80000000#len(balrog)

g = balrog['FLUX_NOISELESS_G'][:n_try]
r = balrog['FLUX_NOISELESS_R'][:n_try]
i = balrog['FLUX_NOISELESS_I'][:n_try]
z = balrog['FLUX_NOISELESS_Z'][:n_try]
y = balrog['FLUX_NOISELESS_Y'][:n_try]
size = balrog['HALFLIGHTRADIUS_0'][:n_try]
index = balrog['BALROG_INDEX'][:n_try]

del balrog

sys.stderr.write('deleted balrog table\n')

vec = np.array([g,r,i,z,y,size])

mag_vec = np.prod(zip(vec, [mu]*(len(vec)-1) + [np.sqrt(mu)]), axis=1)
data = np.array( zip(*[(fi-np.mean(mfi))/np.std(mfi) for fi,mfi in zip(vec,mag_vec)]) )

del vec

st = time.time()
tree = ckdtree.cKDTree(data, balanced_tree=False)
en = time.time()

del data

sys.stderr.write('tree build time ({} elements): {}s\n'.format(n_try, en-st))

mag_data = np.array( zip(*((mfi-np.mean(mfi))/np.std(mfi) for mfi in mag_vec)) )

sys.stderr.write('before query\n')
q_start = time.time()
dists, indices = [], []
for md in mag_data:
    try:
        dist,ind = tree.query(md)
        dists.append(dist)
        indices.append(ind)
    except MemoryError:
        sys.stderr.write('memory error\n')
        printtime()
        exit()
q_end = time.time()
sys.stderr.write('\nquery time: {}s\n'.format(q_end-q_start))

del mag_vec, tree

#print "first element of magnified data: {}, matched vector id: {}".format(mag_data[0], indices[0])

del mag_data

t = Table()
t['BALROG_INDEX'] = index

indices = [int(ii) for ii in indices]

t['BALROG_INDEX_MAG'+str(mu)] = t['BALROG_INDEX'][indices]
t['FLUX_NOISELESS_G_MAG'+str(mu)] = g[indices]
t['FLUX_NOISELESS_R_MAG'+str(mu)] = r[indices]
t['FLUX_NOISELESS_I_MAG'+str(mu)] = i[indices]
t['FLUX_NOISELESS_Z_MAG'+str(mu)] = z[indices]
t['FLUX_NOISELESS_Y_MAG'+str(mu)] = y[indices]
t['HALFLIGHTRADIUS_0_MAG'+str(mu)] = size[indices]
#t['D_MAG'+str(mu)] = dists
#t['I_MAG'+str(mu)] = indices

print "output table length: {}".format(len(t))

del g,r,i,z,y,index,size,dists,indices

if os.path.exists(output):
    spl  = os.path.splitext(output)
    output = "{}-new{}".format(spl[0],spl[1])
    print "output file exists: output will be written to {}".format(output)
    
t.write(output)
print "output written to {}".format(output)

printtime()
