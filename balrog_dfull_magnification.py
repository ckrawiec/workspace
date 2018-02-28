import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import glob
from astropy.io import fits
from scipy.spatial import ckdtree
from scipy.interpolate import griddata

#remain in linear limit of dn/dmu
#take sizes into account as well
#fluxes are different magnitudes in different filters
#size in different units than flux, log size
#can use magnitudes, but negative flux can change mag significantly
#lower limit in g magnitude, S/N above limit, 2-3
#can do logs, cut on dfull low flux, same limit on noiseless

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


def createdata(dfull, sdfull, sizes=True, logs=True):
    if logs:
        df_vec = [np.log10(dfull['FLUX_AUTO_'+f]) for f in filters]
        #magnified fluxes
        new_df_vec = df_vec + np.log10(mu)
    else:
        df_vec = [dfull['FLUX_AUTO_'+f] for f in filters]
        new_df_vec = df_vec * mu

    dfull_hlr = gethlr(dfull, sdfull)
    if sizes:
        df_vec = np.vstack([df_vec, dfull_hlr])
        
        #magnified sizes
        new_dfull_hlr = dfull_hlr * np.sqrt(mu)
        new_df_vec = np.vstack([new_df_vec, new_dfull_hlr])
        
    df_data = np.array( zip(*df_vec) )[dfull_hlr>0.]
    
    new_df_data = np.array( zip(*new_df_vec) )[dfull_hlr>0.]
    return df_data, new_df_data

def createbalrog(brog, sizes=True, logs=True):
    if logs:
        br_vec = [np.log10(brog['FLUX_NOISELESS_'+f]) for f in filters]
    else:
        br_vec = [brog['FLUX_NOISELESS_'+f] for f in filters]
    if sizes:
        br_vec = np.vstack([br_vec, brog['HALFLIGHTRADIUS_0']])
    br_data = np.array( zip(*br_vec) )
    return br_data
    

def gethlr(dfull, sdfull):
    #linear fit using just flux radii
    #(dfull['FLUX_RADIUS_I'] - 1.7) / 2.3
    
    #use average star flux_radius_i around galaxies
    #& galaxy flux_radius from balrog
    to_grid = fits.open('/Users/Christina/DES/data/balrog/sva1/balrog_tab01_avg_star_fluxradiusi_0.1deg.fits')[1].data

    #make tree of dfull stars
    sys.stderr.write("stars in dfull: {}\n".format(len(sdfull)))
    sys.stderr.write("galaxies in dfull: {}\n".format(len(dfull)))
    
    sys.stderr.write('creating hlr tree...\n')
    star_tree = ckdtree.cKDTree(zip(sdfull['RA'], sdfull['DEC']))
    
    gal_pos = zip(dfull['RA'], dfull['DEC'])
    #360 arcsec
    close = star_tree.query_ball_point(gal_pos, r=0.1)
    sys.stderr.write('    calculating average star radii...\n')
    
    start = time.time()
    dfull_avg_star_fr = np.array([np.median(sdfull['FLUX_RADIUS_I'][c]) for c in close])
    end = time.time()
    
    dfull_hlr = griddata(zip(to_grid['flux_radius_i'], to_grid['avg_flux_radius_i']), to_grid['hlr'],
                        zip(dfull['FLUX_RADIUS_I'], dfull_avg_star_fr))
    sys.stderr.write('    done\n')
    return dfull_hlr

def getslope(table, column, mask):
    h = np.histogram(table[column][mask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1,len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    slope = 2.5 * b

    return slope

def plotslopes(m, a, t, r, s):
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

def main():
    dfull = fits.open('/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits')[1].data

    #cut out low fluxes
    flux_mask = (dfull['FLUX_AUTO_G'] > cut) & (dfull['FLUX_AUTO_R'] > cut) & (dfull['FLUX_AUTO_I'] > cut) & (dfull['FLUX_AUTO_Z'] > cut) & (dfull['FLUX_RADIUS_I'] > 0.)

    #separate galaxies & stars
    dfull_gals = dfull['MODEST_CLASS']==1
    dfull_stars = dfull['MODEST_CLASS']==0
    sdfull = dfull[dfull_stars]
    dfull = dfull[dfull_gals & flux_mask]

    df_data, new_df_data = createdata(dfull, sdfull, sizes=True, logs=True)

    #indices of dfull data
    ids = list(range(len(df_data)))

    #cycle over all balrog tables #should use zpcorr ones??????
    all_d, all_new_d = [], []
    bi, new_bi =[], []
    tabnums = []
    #up to data/ is 25
    for tab in glob.glob('/Users/Christina/DES/data/balrog_sva1_tab*_TRUTH_zp_corr_fluxes.fits'):
        tabnum = tab[41:43]
        tabnums.append(tabnum)
        sys.stderr.write('working on balrog table {}...\n'.format(tabnum))

        #open tables, get galaxies
        brog = fits.open(tab)[1].data
        
        br_flux_mask = (brog['FLUX_NOISELESS_G'] > 0.) & (brog['FLUX_NOISELESS_R'] > 0.) & (brog['FLUX_NOISELESS_I'] > 0.) & (brog['FLUX_NOISELESS_Z'] > 0.) & (brog['HALFLIGHTRADIUS_0'] > 0.)

        #choose just galaxies #???
        brog_gals = brog['OBJTYPE']==1
        brog = brog[brog_gals & br_flux_mask]
        br_data = createbalrog(brog, sizes=True, logs=True)
        
        #create balrog tree
        sys.stderr.write('    creating tree...\n')
        br_tree = ckdtree.cKDTree(br_data)
    
        #query tree for fluxes
        sys.stderr.write('    querying original...\n')
        orig_d, orig_id = br_tree.query(df_data)
        all_d.append(orig_d)
        
        #query tree for magnified fluxes
        sys.stderr.write('    querying magnified...\n')
        new_d, new_id = br_tree.query(new_df_data)
        all_new_d.append(new_d)
        
        #save d's, need to choose closest one later
        #find balrog ids
        orig_index = brog['BALROG_INDEX'][orig_id]
        new_index = brog['BALROG_INDEX'][new_id]
        bi.append(orig_index)
        new_bi.append(new_index)

    best_tabs = np.argmin(zip(*all_d), axis=1)
    new_best_tabs = np.argmin(zip(*all_new_d), axis=1)
    orig_olap, new_olap = 0, 0
    for itab in range(len(tabnums)):
        #open sim catalog
        sim = fits.open('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM.fits'.format(tabnums[itab]))[1].data

        #focus on this table
        this_set = np.where(best_tabs==itab)
        this_new_set = np.where(new_best_tabs==itab)
        
        #count objects that are found (detected) in sim catalog
        sys.stderr.write('working on sim table {}...\n'.format(tabnums[itab]))
        sys.stderr.write('    checking for detections...\n')

        all_d[itab][this_set]
        #detected objects whose truth fluxes match dfull
        orig_olap += len(set(sim['BALROG_INDEX']).intersection(bi[itab][this_set]))
        
        #detected objects whose truth fluxes match magnified dfull
        new_olap += len(set(sim['BALROG_INDEX']).intersection(new_bi[itab][this_new_set]))
        
    #report information
    #print "Matches with original flux: {}".format(len(orig_d))
    #print "Matches with flux*{}: {}".format(mu, len(new_d))
    
    print "Detected original matches: {}".format(orig_olap)
    print "Detected magnified matches: {}".format(new_olap)
    
    #totally doing this wrong => check sim magnitudes for window matching
    """
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
        t.append(getslope(dfull, 'MAG_AUTO_I', mask))
    
        #balrog truth
        r.append(getslope(brog, 'MAG_I', brmask))
    
        #balrog sim
        s.append(getslope(sim, 'MAG_AUTO_I', smask))

    plotslopes(m, a, t, r, s)
    """
if __name__=="__main__":
    main()    
