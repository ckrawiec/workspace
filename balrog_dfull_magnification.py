import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.io import fits
from scipy.spatial import ckdtree

#remain in linear limit of dn/dmu
#take sizes into account as well
#fluxes are different magnitudes in different filters
#size in different units than flux, log size
#can use magnitudes, but negative flux can change mag significantly
#lower limit in g magnitude, S/N above limit, 2-3

filters = 'GRIZ'
mu = 1.1

cut = 40.

#try a few different bins, see that shape of slope follows
mag_limits = [[20., 21.],
              [20.5, 21.5],
              [21., 22.],
              [21.5, 22.5],
              [22., 23.],
              [22.5, 23.5],
              [23., 24.],
              [23.5, 24.5],
              [24., 25.]]

#open tables, get galaxies
brog = fits.open('/Users/Christina/DES/data/balrog_sva1_tab01_TRUTH_fluxes.fits')[1].data
dfull = fits.open('/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits')[1].data

#cut out low fluxes
flux_mask = (dfull['FLUX_AUTO_G'] > cut) & (dfull['FLUX_AUTO_R'] > cut) & (dfull['FLUX_AUTO_I'] > cut) & (dfull['FLUX_AUTO_Z'] > cut)

#choose just galaxies
brog_gals = brog['OBJTYPE']==1
dfull_gals = dfull['MODEST_CLASS']==1
brog = brog[brog_gals]
dfull = dfull[dfull_gals & flux_mask]

dfull_hlr = (dfull['FLUX_RADIUS_I'] - 1.7) / 2.3

df_vec = [np.log10(dfull['FLUX_AUTO_'+f]) for f in filters]
#df_vec.append(list(dfull_hlr))
br_vec = [np.log10(brog['FLUX_NOISELESS_'+f]) for f in filters]
#br_vec.append(list(brog['HALFLIGHTRADIUS_0']))

#can do logs, cut on dfull low flux, same limit on noiseless
df_data = np.array( zip(*df_vec) )
#np.array( zip(*[dfull['FLUX_AUTO_'+f]/np.median(dfull['FLUX_AUTO_'+f]) for f in filters]) )
br_data = np.array( zip(*br_vec) )
#np.array( zip(*[brog['FLUX_NOISELESS_'+f]/np.median(dfull['FLUX_AUTO_'+f]) for f in filters]) )

#magnified data
new_df_data = df_data + np.array([np.log10(mu)]*4)# + [0])
#new_df_data = df_data * np.array([1]*4 + [np.sqrt(mu)])

#indices of dfull data
ids = list(range(len(df_data)))

#create balrog tree
sys.stderr.write('creating tree...\n')
br_tree = ckdtree.cKDTree(br_data)

#query tree for fluxes
sys.stderr.write('querying...\n')
orig_d, orig_id = br_tree.query(df_data)

#query tree for magnified fluxes
sys.stderr.write('querying...\n')
new_d, new_id = br_tree.query(new_df_data)

#report information
print "Matches with original flux: {}".format(len(orig_d))
print "Matches with flux*{}: {}".format(mu, len(new_d))

print "All matches the same? : \n", np.all(orig_id==new_id)

#find balrog ids
orig_index = brog['BALROG_INDEX'][orig_id]
new_index = brog['BALROG_INDEX'][new_id]


#open sim catalog
sim = fits.open('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab01_SIM.fits')[1].data
#balrog_sva1_auto_tab01_SIM_TRUTH_zp_corr_fluxes.fits')[1].data

#count objects that are found (detected) in sim catalog
sys.stderr.write('checking for detections...\n')

#detected objects whose truth fluxes match dfull
orig_olap = set(sim['BALROG_INDEX']).intersection(orig_index)
#detected objects whose truth fluxes match magnified dfull
new_olap = set(sim['BALROG_INDEX']).intersection(new_index)

print "Detected original matches: {}/{}".format(len(orig_olap), len(orig_index))
print "Detected magnified matches: {}/{}".format(len(new_olap), len(new_index))

#totally doing this wrong => check sim magnitudes for window matching

#check in magnitude range
m, a, t, r, s = [], [], [], [], []
for mag_min, mag_max in mag_limits:
    mask = (dfull['MAG_AUTO_I'] > mag_min) & (dfull['MAG_AUTO_I'] < mag_max)
    brmask = (brog['MAG_I'] > mag_min) & (brog['MAG_I'] < mag_max)
    smask = (sim['MAG_AUTO_I'] > mag_min) & (sim['MAG_AUTO_I'] < mag_max)
    
    orig_olap_in_mag = len(set(sim['BALROG_INDEX'][smask]).intersection(orig_index))
    new_olap_in_mag = len(set(sim['BALROG_INDEX'][smask]).intersection(new_index))

    slope = (new_olap_in_mag - float(orig_olap_in_mag)) / mu
    m.append(np.mean([mag_min, mag_max]))
    try:
        a.append((slope + orig_olap_in_mag) / orig_olap_in_mag)
    except ZeroDivisionError:
        a.append(np.nan)
    print "Detected originals in i_mag range {}-{}: {}".format(mag_min, mag_max, orig_olap_in_mag)
    print "Detected magnified in i_mag range {}-{}: {}\n".format(mag_min, mag_max, new_olap_in_mag)

    #dfull
    h = np.histogram(dfull['MAG_AUTO_I'][mask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1,len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    t.append(2.5 * b)

    #balrog truth
    h = np.histogram(brog['MAG_I'][brmask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1, len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    r.append(2.5 * b)

    #balrog sim
    h = np.histogram(sim['MAG_AUTO_I'][smask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1,len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    s.append(2.5 * b)
    

#plot measured slopes
plt.plot(m, a, label='mine')
#plot actual slopes
plt.plot(m, t, label='dfull auto')
#plot balrog slopes
plt.plot(m, r, label='balrog noiseless')
plt.plot(m, s, label='balrog sim')

plt.ylabel('2.5 * d(log n(m))/dm')
plt.xlabel('i magnitude')
plt.legend(loc='best')

plt.ylim(-0.5, 1.2)
plt.savefig('/Users/Christina/DES/magnification/balrog_dfull_magnification')
