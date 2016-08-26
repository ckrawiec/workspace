#a few different ways to add noise
# - get noise from neighbors for a given ra/dec
# - get limiting mag from map, use mag/err relation
import pandas
import itertools
import time
import os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import ascii, fits
from multiprocessing import Pool
from zprobability import p, ptree

data_file = '/Users/Christina/DES/data/sva1_gold_detmodel_gals.fits'
base_data_type = 'mag_detmodel_'
base_err_type = 'magerr_detmodel_'

filters = ['g','r','i','z','Y']

sva1_cosmos_file = '{}/DES/data/sva1_coadd_cosmos.fits'.format(home_dir)
cosmos_file = '{}/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'.format(home_dir)

#indices for sva1_gold_detmodel_gals.fits                                                                                    
gold_cosmos15_indices_file = '{}/DES/magnification/lbgselect/gold_cosmos15_indices.txt'.format(home_dir)
gold_no_cosmos_indices_file = '{}/DES/magnification/lbgselect/gold_no_cosmos_indices.txt'.format(home_dir)
cosmos15_indices_file = '{}/DES/magnification/lbgselect/cosmos15_indices.txt'.format(home_dir)

ptype = 'full'

def pwrapper(args):
    if ptype=='tree':
        return ptree(*(args+[knear=k_near]))
    elif ptype=='full':
        return p(*args)
    else:
        raise ValueError('Choose ptype=\'full\' or \'tree\'')

def main():
    #open data catalog
    data = fits.open(data_file)[1].data
    
    ra, dec = data['ra'], data['dec']

    #healpixify errors (mean or median?)    
    theta = (90.0 - dec)*np.pi/180.
    phi = ra*np.pi/180.
    
    nside = 512
    pix = hp.ang2pix(nside, theta, phi, nest=True)
    hpmap = np.zeros(hp.nside2npix(nside), dtype=np.double)
    
    noisemaps = {}
    for filter in filters:
        df = pandas.DataFrame( {'data': data[base_data_type+filter],
                                'err': data[base_err_type+filter],
                                'pixel': pix} )
    
        gpix = df.groupby(['pixel'])
        meanerr = gpix['err'].mean()
        noisemaps[filter] = meanerr
        
    gold_cosmos15 = np.loadtxt(gold_cosmos15_indices_file, dtype=int)
    gold_no_cosmos = np.loadtxt(gold_no_cosmos_indices_file, dtype=int)
    cosmos15_indices = np.loadtxt(cosmos15_indices_file, dtype=int)

    #COSMOS2015 photo-z                                    
    #z_cosmos = 9.99 --> X-ray object, z_cosmos = 0 --> star
    z_cosmos = cosmos15['photoz'][cosmos15_indices]
    z0mask = (z_cosmos > 0) & (z_cosmos < 9.9)
    z3mask = (z_cosmos >= 3.) & (z_cosmos < 9.9)
    z4mask = (z_cosmos >= 4.) & (z_cosmos < 9.9)

    Ntot = len(z_cosmos[z0mask & ~z4mask]) + len(z_cosmos[z4mask])

    #cosmos fluxes and errors from sva1 gold                                                                                   
    def maketable(mask=None, cosmos=False, filters=['g','r','i','z','Y']):
        table = {}
        table[f] = sva1_gold[data_type+'_detmodel_'+f][gold_cosmos15][mask]
        table[f+'err'] = sva1_gold[data_type+'err_detmodel_'+f][gold_cosmos15][mask]
        return table




#open cosmos catalog
#separate into groups
#at each RA/DEC pixel
# for each groupi
#  for each groupj
#   compute P(groupj|groupi+noise) 

#just a subsample of whole group? Exclude that galaxy/ies from the group?

#save results to table
