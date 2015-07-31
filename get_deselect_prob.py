import numpy as np
from astropy.io import fits
import glob

moment_files = glob.glob(sys.argv[1])

deselect_prob_list = []

for moment_file in moment_files:
    f = fits.open(moment_file)

    d = f[1].data

    nuniq = d['NUNIQUE']
    pqr = d['PQR']

    i = 0
    while i < len(nuniq):
        if nuniq[i]==-1:
            t = pqr[i]
            p = t[0]
            q = np.array([t[1],t[2]])
            r = np.matrix([t[3],t[4]],
                          [t[4],t[5]])
            break
        i+=1

    deselect_prob_list.append(p + np.inner(q,g) + 0.5*np.inner(g,np.inner(r,g)))

deselect_prob = np.mean(deselect_prob_list)
select_prob = 1-deselect_prob

u = np.std(deselect_prob)/sqrt(len(moment_files))

print "Probability of deselection: ", deselect_prob
print "Probability of selection: ", select_prob
print "Uncertainty: ", u




