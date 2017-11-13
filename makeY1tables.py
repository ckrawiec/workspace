import numpy as np
import healpy as hp
import sys
import glob
import os
import esutil
import matplotlib.pyplot as plt
from astropy.table import Table, Column, vstack, join

home_dir = '/Users/Christina/'#'/home/ckrawiec/'

y1_gold_z3_4bins_highP = '/Users/Christina/DES/magnification/lbgselect/y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_4bins_full_gauss_000001_Phigh0.99_positions.csv'

y1a1_gold_tables = glob.glob(home_dir+'DES/data/y1a1_gold_flux_auto_griz_*')
y1a1_gold = home_dir+'DES/data/y1a1_gold_flux_auto_griz.fits'
y1a1_d04_tables = glob.glob(home_dir+'DES/data/y1a1_gold_d04_0000*fits')
y1a1_d10_tables = glob.glob(home_dir+'DES/data/y1a1_gold_d10_0000*fits')

y1d04dfullzp = '/Users/Christina/DES/magnification/lbgselect/y1a1_gold_cosmosd04_cosmosdfull_matched_zminchi2_auto_griz_full_z3_4bins_catmatch_both.fits'
y1_zp_tables = glob.glob(home_dir+'DES/magnification/lbgselect/zproboutput/y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_4bins_full_gauss_000001*fits')
y1_zp_table = home_dir+'DES/magnification/lbgselect/zproboutput/y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_4bins_full_gauss_000001.fits'
y1_zp_match = '/Users/Christina/DES/magnification/lbgselect/y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_4bins_full_gauss_000001_catmatch.fits'

y1_dfull = home_dir+'DES/data/y1a1_gold_dfull.fits'
y1_d04 = home_dir+'DES/data/y1a1_gold_d04.fits'
y1_d10 = home_dir+'DES/data/y1a1_gold_d10.fits'

y1_dfull_cosmos_masked = home_dir+'DES/data/y1a1_gold_dfull_masked_cosmos.fits'
y1_dfull_masked = home_dir+'DES/data/y1a1_gold_dfull_masked.fits'

y1_dfull_match_file = home_dir+'DES/data/match_y1a1_gold_dfull_cosmos_1arcsec'
y1_dfull_masked_cosmos_match_file = home_dir+'DES/data/match_y1a1_gold_dfull_masked_cosmos_1arcsec'
y1_d04_match_file = home_dir+'DES/data/match_y1a1_gold_d04_cosmos_1arcsec'
d04_dfull_cosmos_match_file = home_dir+'DES/data/match_y1a1_gold_d04_dfull_cosmos_1arcsec'
d04_dfull_cosmos_masked_match_file = home_dir+'DES/data/match_y1a1_gold_d04_dfull_masked_cosmos_1arcsec'
d10_dfull_cosmos_match_file = home_dir+'DES/data/match_y1a1_gold_d10_dfull_cosmos_1arcsec'
d04_dfull_match_file = home_dir+'DES/data/match_y1a1_gold_d04_dfull_1arcsec'
d10_dfull_match_file = home_dir+'DES/data/match_y1a1_gold_d10_dfull_1arcsec'

y1_d04_cosmos = home_dir+'DES/data/y1a1_gold_d04_cosmos.fits'
y1_dfull_cosmos = home_dir+'DES/data/y1a1_gold_dfull_cosmos.fits'
y1_d10_cosmos = home_dir+'DES/data/y1a1_gold_d10_cosmos.fits'

d04_dfull_cosmos = home_dir+ 'DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits'
d04_dfull_cosmos_masked = home_dir + 'DES/data/y1a1_gold_d04_dfull_masked_cosmos_matched.fits'
d10_dfull_cosmos = home_dir+'DES/data/y1a1_gold_d10_dfull_cosmos_matched.fits'
d04_dfull = home_dir+'DES/data/y1a1_gold_d04_dfull_matched.fits'

cosmos_file = home_dir+'COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'

balrog_glob = home_dir + 'DES/data/balrog_y1a1_truth_index_z_*'
balrog_output = home_dir + 'DES/data/balrog_y1a1_truth_index_z.fits'

ngmix_dfull = home_dir + 'DES/data/DES0215-0458-y1a1-DFULL-mof-crucial-001.fits'
ngmix_dfull_match = os.path.splitext(ngmix_dfull)[0]+'_match_cosmos_1arcsec'

junk_match_file = home_dir+'junk_match'

def maketables():
#    makebalrog()
#    chooseobjtype1(balrog_output)
#    matchwithcosmos(ngmix_dfull, ngmix_dfull_match)
#    stackwrite(y1a1_gold_tables, y1a1_gold)
#    stackwrite(y1_zp_tables, y1_zp_table)
#    mask(y1_zp_match)
#    flagmask(y1_zp_match)
#    flagmask(y1d04dfullzp, col='FLAGS_GOLD_d04')
#    flagmask(d04_dfull_cosmos, col='FLAGS_GOLD_d04', tag='d04flagmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagmasked.fits',
#             col='FLAGS_GOLD_dfull', tag='dfullflagmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits',
#             col='FLAGS_BADREGION_d04', tag='d04flagBRmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBRmasked.fits',
#             col='FLAGS_BADREGION_dfull', tag='dfullflagBRmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBRmasked_dfullflagBRmasked.fits',
#             col='FLAGS_GOLD_d04', tag='d04flagmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBRmasked_dfullflagBRmasked_d04flagmasked.fits',
#             col='FLAGS_GOLD_dfull', tag='dfullflagmasked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits',
#             col='FLAGS_BADREGION_d04', val=64, tag='d04flagBR64masked')
#    flagmask('/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04flagBR64masked.fits',
#             col='FLAGS_BADREGION_dfull', val=64, tag='dfullflagBR64masked')
    flagmask('/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits',
             col='FLAGS_BADREGION', tag='flagBR64masked')
#    matchwithcosmos(y1_d04, y1_d04_match_file)
#    chooseobjtype1(y1_d04)
#    chooseobjtype1(y1_d04_cosmos)
#    mask(y1_dfull, y1_dfull_masked)
#    mask(d04_dfull_cosmos, ra_col='RA_d04', dec_col='DEC_d04')
#    mask(home_dir+ 'DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits',
#         ra_col='RA_d04', dec_col='DEC_d04',
#         tag='d04masked', mask_map='y1d04')
#    mask(home_dir+ 'DES/data/y1a1_gold_d04_dfull_cosmos_matched_d04masked.fits',
#         ra_col='RA_dfull', dec_col='DEC_dfull',
#         tag='dfullmasked', mask_map='y1dfull')
#    mask(y1_gold_z3_4bins_highP)
#    mask(y1a1_gold)
#    stackwrite(y1a1_d10_tables, y1_d10)
#    match(y1_d10, y1_dfull_cosmos, 'd10', 'dfull', d10_dfull_cosmos_match_file, d10_dfull_cosmos)
#    match(y1_d04, y1_dfull_cosmos_masked, 'd04', 'dfull', d04_dfull_cosmos_masked_match_file, d04_dfull_cosmos_masked)
#    match(y1d04dfullzp, y1_dfull_cosmos, 'd04', 'dfull', junk_match_file, os.path.splitext(y1d04dfullzp)[0]+'_both.fits')
#    matchwithcosmos(y1_dfull_masked, y1_dfull_masked_cosmos_match_file)
#    separatecosmosstars(y1_dfull_cosmos_masked)
#    chunkandwrite(d04_dfull_cosmos, 2)


def separatecosmosstars(table):
    #separate table into stars vs no stars by cosmos type
    sys.stderr.write('identifying stars using COSMOS...\n')
    tab = Table.read(table)
    stars = tab['TYPE']==1
    star_tab = tab[stars]
    star_tab.write(os.path.splitext(table)[0]+'_stars.fits')

    no_star_tab = tab[~stars]
    no_star_tab.write(os.path.splitext(table)[0]+'_no_stars.fits')
    sys.stderr.write('done\n\n')

def chunkandwrite(table, nchunks):
    sys.stderr.write('dividing table...\n')
    tab = Table.read(table)
    inds = np.arange(len(tab))
    chunks = np.array_split(inds, nchunks)
    ichunk = 1
    for chunk in chunks:
        sys.stderr.write('    working on chunk {}\n'.format(ichunk))
        newtab = tab[chunk]
        newtab.write(os.path.splitext(table)[0]+'_chunk{}-{}.fits'.format(str(ichunk).zfill(2), nchunks))
        ichunk+=1
    sys.stderr.write('done\n\n')

def stackwrite(flist, output_file):
    tlist = [Table.read(table) for table in flist]
    new_table = vstack(tlist)
    del tlist
    new_table.write(output_file)

def flagmask(table, col='FLAGS_GOLD', val=None, tag='flagmasked'):
    tab = Table.read(table)
    if val is not None:
        new_tab = tab[(tab[col]==0) | ((tab[col]-val)==0)]
    else:
        new_tab = tab[tab[col]==0]
    output_file = os.path.splitext(table)[0]+'_{}.fits'.format(tag)
    sys.stderr.write('writing masked file to {}...\n'.format(output_file))
    new_tab.write(output_file)
    
def mask(table, ra_col='RA', dec_col='DEC', tag='masked', mask_map=''):
    if (mask_map=='y1dfull'):# & (table==y1_dfull_cosmos) or (table==y1_dfull):
        sys.stderr.write("using map 'y1dfull'\n")
        footmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.2_dfull_footprint_4096.fits.gz')
        badmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.3_dfull_badmask_4096.fits.gz',
                            dtype=np.int32)
    elif (mask_map=='y1gold'):# or (table==y1a1_gold) or (table==y1_gold_z3_4bins_highP) or (y1_zp_match):
        sys.stderr.write("using map 'y1gold'\n")
        footmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.2_wide_footprint_4096.fits.gz')
        badmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.3_wide_badmask_4096.fits.gz',
                            dtype=np.int32)
    elif (mask_map=='y1d04'):#(table==d04_dfull_cosmos):
        sys.stderr.write("using map 'y1d04'\n")
        footmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.2_d04_footprint_4096.fits.gz')
        badmask=hp.read_map(home_dir+'DES/data/y1a1_gold_1.0.3_d04_badmask_4096.fits.gz',
                            dtype=np.int32)
    else:
        return False

    tab = Table.read(table)
    print "Table length before masking: ", len(tab)
    nside = hp.npix2nside(footmask.size)
    print "NSIDE = ", nside

    theta = (90.0 - tab[dec_col])*np.pi/180.
    phi = tab[ra_col]*np.pi/180.
    pix = hp.ang2pix(nside, theta, phi)
    ipring, = np.where((footmask[pix] >= 1) & (badmask[pix] == 0))

    new_table = tab[ipring]
    print "Table length after masking: ", len(new_table)

    del tab
    sys.stderr.write("deleted old table, writing to file...\n")

    output_file = os.path.splitext(table)[0]+'_{}.fits'.format(tag)
    new_table.write(output_file)
    

def match(file1, file2, name1, name2, match_file, outfile):
    sys.stderr.write('matching {} and {}...\n'.format(file1, file2))
    tab1 = Table.read(file1)
    tab2 = Table.read(file2)

    sys.stderr.write('    tables read in\n')

    h = esutil.htm.HTM(10)
    h.match(tab1['RA'], tab1['DEC'],
            tab2['RA'], tab2['DEC'],
            radius=1./3600,
            file=match_file)

    sys.stderr.write('    matched\n')

    m = h.read(match_file)

    m1, m2, merr = np.array(zip(*m))

    sys.stderr.write('    {} matches\n'.format(len(m1)))
    sys.stderr.write('    No duplicates in matched list?: {}\n'.format(len(m1)==len(set(m1))))

    m1_int = [int(g) for g in m1]
    m2_int = [int(c) for c in m2]

    new_tab = Table()
    for colname in tab1.colnames:
        new_tab.add_column(Column(data=tab1[colname][m1_int], name=colname+'_'+name1))
    for colname in tab2.colnames:
        new_tab.add_column(Column(data=tab2[colname][m2_int], name=colname+'_'+name2))

    new_tab.add_column(Column(name='match_err', data=merr))
    new_tab.write(outfile)
    sys.stderr.write('done - matched table written to {}\n\n'.format(outfile))

def matchwithcosmos(tab_file, match_file, nocosmos=False):
    print "Matching {} with {} within 1 arcsec...".format(tab_file, cosmos_file)

    tab = Table.read(tab_file)

    if 'mof' in tab_file:
        y1tab = Table.read(y1_dfull)
        ntab.rename_column('id', 'COADD_OBJECTS_ID')
        tab = join(tab, y1tab, keys='COADD_OBJECTS_ID')

        print ntab['COADD_OBJECTS_ID'][:10]
        print len(tab)
        print y1tab['COADD_OBJECTS_ID'][:10]
        plt.scatter(tab['RA'], tab['DEC'])
        plt.show()

    cosmos = Table.read(cosmos_file)

    h = esutil.htm.HTM(10)
    h.match(tab['RA'], tab['DEC'],
            cosmos['ALPHA_J2000'], cosmos['DELTA_J2000'],
            radius=1./3600,
            file=match_file)

    m = h.read(match_file)

    data_m, cosmos_m, merr = np.array(zip(*m))
               
    print " {} COSMOS/Y1 matches".format(len(cosmos_m))
    print " No duplicates in matched list?:{}".format(len(data_m)==len(set(data_m)))

    if nocosmos==True:
        no_cosmos = list(set(range(0,len(tab)))-set(data_m))
        new_data = tab[no_cosmos]
        print "made new table without cosmos matches"
        new_data.write(os.path.splitext(tab_file)[0]+'_no_cosmos.fits')

        del new_data

        print ' wrote new table'

    data_m_int = [int(g) for g in data_m]
    cosmos_m_int = [int(c) for c in cosmos_m]

    new_cosmos = tab[data_m_int] 
    del tab

    print ' deleted old table'

    new_cosmos.add_column(Column(name='photoz', data=cosmos['PHOTOZ'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='ZMINCHI2', data=cosmos['ZMINCHI2'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='TYPE', data=cosmos['TYPE'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='NUMBER', data=cosmos['NUMBER'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='ALPHA_J2000', data=cosmos['ALPHA_J2000'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='DELTA_J2000', data=cosmos['ALPHA_J2000'][cosmos_m_int]))
    new_cosmos.add_column(Column(name='match_err', data=merr))

    print ' made new cosmos table'

    new_cosmos.write(os.path.splitext(tab_file)[0]+'_cosmos.fits')


def makebalrog():
    balrog_tables = glob.glob(balrog_glob)
    stackwrite(balrog_tables, output)

def chooseobjtype1(output):
    t = Table.read(output)
    tnew = t[t['OBJTYPE']==1]
    del t
    tnew.write(os.path.splitext(output)[0]+'_objtype1.fits')
    
maketables()

#*********Y1A1 GOLD DFULL********

#tlist = []
#for table in y1a1_gold_tables:
#    tlist.append(Table.read(table))
#new_Y1 = vstack(tlist)
#del tlist
#new_Y1.write('/home/ckrawiec/DES/data/y1a1_gold_flux_detmodel_MC1.fits')
#del new_Y1
#y1a1_gold_dfull_tables = glob.glob('/home/ckrawiec/DES/data/y1a1_gold_dfull_*')
#tlist = []
#for table in y1a1_gold_dfull_tables:
#    tlist.append(Table.read(table))
#new_Y1_dfull = vstack(tlist)
#del tlist
#new_Y1_dfull.write('/home/ckrawiec/DES/data/y1a1_gold_dfull.fits')



