import matplotlib.pyplot as plt
from astropy.table import Table, join
from astropy.io import fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


inputfile  = '/Users/Christina/Galaxia/output/DESSV/SV390deg.ebf'
outputfile = '/Users/Christina/Galaxia/output/DESSV/SV390deg_desflux.fits'

#SDSS->DES transformations from Douglas Tucker
#The u, g, r transformations apply for stars with 0.2 <= (g-r)_sdssdr13 < 1.2.
#The i, z, Y transformations apply for stars with 0.0 <= (i-z)_sdssdr13 < 0.8.

def des(g,r,i,z,u=None):
    
    g_des = g + 0.001 - 0.075*(g-r)                   #(RMS=0.021 mag per star)
    r_des = r - 0.009 - 0.069*(g-r)                   #(RMS=0.021 mag per star)
    i_des = i + 0.014 - 0.214*(i-z) - 0.096*(i-z)**2  #(RMS=0.023 mag per star)
    z_des = z + 0.022 - 0.068*(i-z)                   #(RMS=0.025 mag per star)
    Y_des = z + 0.045 - 0.306*(i-z)                   #(RMS=0.030 mag per star)
    if u:
        u_des = u - 0.479 + 0.466*(g-r) - 0.350*(g-r)**2  #(RMS=0.055 mag per star)
        return g_des,r_des,i_des,z_des,u_des
    return g_des,r_des,i_des,z_des

def galaxia_to_fits(galtab, filename):
    #convert to apparent mags, then to DES, get fluxes, and save new table
    for filt in 'griz':
        galtab['app_sdss_'+filt] = galtab['sdss_'+filt] + 5*np.log10(100*galtab['rad'])
    galtab['app_des_g'], galtab['app_des_r'], galtab['app_des_i'], galtab['app_des_z'] = des(galtab['app_sdss_g'], 
                                                                                             galtab['app_sdss_r'], 
                                                                                             galtab['app_sdss_i'],
                                                                                             galtab['app_sdss_z'])
    for filt in 'griz':
        galtab['flux_des_'+filt] = 10.**(-0.4*(galtab['app_des_'+filt]-30.))
    
    galtab.write(filename)



def main():
    newgale = ebf.read(inputfile)
    newgale.pop('log')
    newgale.pop('center')
    newgal = Table(newgale)

    galaxia_to_fits(newgal, outputfile)

if __name__=="__main__":
    main()
