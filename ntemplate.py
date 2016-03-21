import numpy as np
from astropy.io import fits
import glob
import sys

moment_files = glob.glob(sys.argv[1])

mins, maxs = [],[]
umins, umaxs = [],[]

for moment_file in moment_files:
    f = fits.open(moment_file)

    d = f[1].data

    used = np.where(d['nunique']!=-1)

    ntemp = np.array(d['ntemplate'])
    mf = np.array([m[0] for m in d['moments']])

    mins.append(np.min(mf))
    maxs.append(np.max(mf))
    umins.append(np.min(mf[used]))
    umaxs.append(np.max(mf[used]))
    

print "minimum M_f: ", np.min(mins)
print "average minimum M_f: ", np.mean(mins)
print "maximum M_f: ", np.max(maxs)
print "average maximum M_f: ", np.mean(maxs)
print "minimum used M_f: ", np.min(umins)
print "average minimum used M_f: ", np.mean(umins)
print "maximum used M_f: ", np.max(umaxs)
print "average maximum M_f: ", np.mean(umaxs)



