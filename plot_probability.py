import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from myutils import match

home_dir = '/home/ckrawiec/'

def main():
    sv = fits.open(home_dir+'DES/data/sva1_gold_detmodel_gals.fits')[1].data
    hi = fits.open(home_dir+'DES/magnification/lbgselect/hi-z_mag_probability.fits')[1].data
    lo = fits.open(home_dir+'DES/magnification/lbgselect/lo-z_mag_probability_test.fits')[1].data
    
    hi.sort(), sv.sort(), lo.sort()
    
    svid, rid = match(sv['coadd_objects_id'], hi['coadd_objects_id'])
    
    ri = sv['mag_detmodel_r'][svid] - sv['mag_detmodel_i'][svid]
    gr = sv['mag_detmodel_g'][svid] - sv['mag_detmodel_r'][svid]
    
    
    
    Phi = hi['hi-z_prob'][rid] / (hi['hi-z_prob'][rid] + lo['lo-z_prob'][rid])
    Plo = lo['lo-z_prob'][rid] / (hi['hi-z_prob'][rid] + lo['lo-z_prob'][rid])
    
    print "P(hi-z)>=0.8: {}/{}".format(len(Phi[Phi>=0.8]), len(Phi))
                
    plt.scatter(Plo, Phi, edgecolor='none')
    plt.ylim(0.8, 1.0)
    plt.xlabel('P(lo-z)')
    plt.ylabel('P(hi-z)')
    plt.savefig(home_dir+'DES/magnification/lbgselect/Phi-vs-Plo')
    plt.close()

if __name__=="__main__":
    main()
                               
