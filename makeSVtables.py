import numpy as np
from astropy.table import Table, Column
import healpy as hp
import esutil
import sys
import os

home_dir = '/Users/Christina/'

sva1_gold = home_dir+'DES/data/sva1_gold_auto.fits'
cosmos_file = home_dir+'COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'

brg_pos = home_dir+'DES/data/balrog_sva1_auto_SIM_J2000.fits'

red = home_dir + 'DES/data/redmagic_sva1_public_v6.3_faint.fits'
redmap = home_dir + 'DES/data/redmagic_v6.2.15_fracdet_nside4096_map.fit.gz'

base_output = 'DES/data/sva1_gold_auto'

sva1_regions = home_dir+'DES/data/sva1_gold_r1.0_goodregions_04_n4096.fits.gz'

sva1_gold_cosmos_match = home_dir+'DES/data/match_sva1_gold_cosmos_gals_1arcsec'

def main():
#    goodregions(brg_pos, sva1_regions, ra_col='ALPHAMODEL_J2000_G', dec_col='DELTAMODEL_J2000_G')
#    matchwithcosmos()
#    goodregions('/Users/Christina/DES/data/sva1_gold_ra_dec.fits', sva1_regions)
#    sptetrim('/Users/Christina/DES/data/sva1_gold_ra_dec_good_regions.fits')
#    goodregions('/Users/Christina/DES/data/redmagic_sva1_public_v6.3_faint.fits', sva1_regions)
#    sptetrim('/Users/Christina/DES/data/redmagic_sva1_public_v6.3_faint_good_regions.fits')
#    goodregions('/Users/Christina/DES/data/balrog_sva1_RADEC.fits', sva1_regions)
#    goodregions('/Users/Christina/DES/data/balrog_sva1_auto_tab01_SIM_TRUTH_zp_corr_fluxes.fits', sva1_regions,
#                ra_col='ALPHAMODEL_J2000_G', dec_col='DELTAMODEL_J2000_G')
#    goodregions('/Users/Christina/DES/data/random_sva1_redmagic.fits', sva1_regions)
#    sptetrim('/Users/Christina/DES/data/random_sva1_redmagic_good_regions.fits')
    goodregions('/home/ckrawiec/DES/data/Galaxia_SV390deg_desflux.fits')
    sptetrim('/home/ckrawiec/DES/data/Galaxia_SV390deg_desflux_good_regions.fits')

def sptetrim(tab_file, ra_col='RA', dec_col='DEC'):
    tab = Table.read(tab_file)
    ra_mask = (tab[ra_col] < 95.) & (tab[ra_col] > 60.)
    dec_mask = (tab[dec_col] < -40.) & (tab[dec_col] > -65.)

    new_tab = tab[ra_mask & dec_mask]

    out_file = os.path.splitext(tab_file)[0]+'_SPT-E.fits'
    new_tab.write(out_file)
    
    sys.stderr.write('wrote new table to {}\n'.format(out_file))

def getweights(tab_file, map):
    tab = Table.read(tab_file)

    hpmap = hp.read_map()

def goodregions(tab_file, region_file, ra_col='RA', dec_col='DEC'):
    tab = Table.read(tab_file)

    hpmap = hp.read_map(region_file, nest=True)
    nside = hp.npix2nside(hpmap.size)

    theta = (90.0 - tab[dec_col])*np.pi/180.
    phi = tab[ra_col]*np.pi/180.
    pix = hp.ang2pix(nside, theta, phi, nest=True)
    good, = np.where(hpmap[pix] == 1)

    sys.stderr.write('good regions found\n')
    
    new_tab = tab[good,]
    out_file = os.path.splitext(tab_file)[0]+'_good_regions.fits'
    new_tab.write(out_file)
    
    sys.stderr.write('wrote new table to {}\n'.format(out_file))
    
def matchwithcosmos(tab_file, match_file, ra_col='RA', dec_col='DEC'):
    cosmos = Table.read(cosmos_file)
    tab = Table.read(tab_file)

    h = esutil.htm.HTM(10)
    h.match(tab[ra_col], tab[dec_col],
            cosmos['ALPHA_J2000'], cosmos['DELTA_J2000'],
            radius=1./3600,
            file=match_file)

    m = h.read(match_file)

    gold_m, cosmos_m, merr = np.array(zip(*m))

    sys.stderr.write('matched\n')

    no_cosmos = list(set(range(0,len(tab)))-set(gold_m))

    sys.stderr.write('No duplicates?: {}\n'.format(len(gold_m)==len(set(gold_m))))

    new_sv = tab[no_cosmos]

    out_file = os.path.splitext(tab_file)[0]+'_no_cosmos.fits'
    new_sv.write(out_file)
    del new_sv

    sys.stderr.write('wrote new table to {}\n'.format(out_file))

    gold_m_int = [int(g) for g in gold_m]
    cosmos_m_int = [int(c) for c in cosmos_m]

    new_cosmos = tab[gold_m_int] 
    del tab
    
    sys.stderr.write('deleted old table\n')

    new_cosmos.add_column(Column(name='NUMBER', data=cosmos['NUMBER'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='photoz', data=cosmos['PHOTOZ'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='zminchi2', data=cosmos['ZMINCHI2'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='TYPE', data=cosmos['TYPE'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='match_err', data=merr))

    cosmos_out_file = os.path.splitext(tab_file)[0]+'_cosmos.fits'
    new_cosmos.write(cosmos_out_file)
    
    sys.stderr.write('wrote new table to {}\n'.format(cosmos_out_file))

if __name__=="__main__":
    main()
