import numpy as np
import glob
import os
#import esutil
import matplotlib.pyplot as plt
from astropy.table import Table, Column, vstack, join

home_dir = '/Users/Christina/'#'/home/ckrawiec/'

y1a1_gold_tables = glob.glob('/home/ckrawiec/DES/data/y1a1_gold_flux_auto_griz_*')
y1a1_gold = '/home/ckrawiec/DES/data/y1a1_gold_flux_auto_griz.fits'
y1a1_d04_tables = glob.glob('/home/ckrawiec/DES/data/y1a1_gold_d04_0000*fits')

y1_dfull = home_dir+'DES/data/y1a1_gold_dfull.fits'
y1_d04 = home_dir+'DES/data/y1a1_gold_d04.fits'

y1_dfull_match_file = home_dir+'DES/data/match_y1a1_gold_dfull_cosmos_1arcsec'
y1_d04_match_file = home_dir+'DES/data/match_y1a1_gold_d04_cosmos_1arcsec'
d04_dfull_match_file = home_dir+'DES/data/match_y1a1_gold_d04_dfull_cosmos_1arcsec'

y1_d04_cosmos = home_dir+'DES/data/y1a1_gold_d04_cosmos.fits'
y1_dfull_cosmos = home_dir+'DES/data/y1a1_gold_dfull_cosmos.fits'

d04_dfull = home_dir+ 'DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits'

cosmos_file = home_dir+'COSMOS/data/COSMOS2015_Laigle+_v1.1.fits'

balrog_glob = home_dir + 'DES/data/balrog_y1a1_truth_index_z_*'
balrog_output = home_dir + 'DES/data/balrog_y1a1_truth_index_z.fits'

ngmix_dfull = home_dir + 'DES/data/DES0215-0458-y1a1-DFULL-mof-crucial-001.fits'
ngmix_dfull_match = os.path.splitext(ngmix_dfull)[0]+'_match_cosmos_1arcsec'

def maketables():
#    makebalrog()
#    chooseobjtype1(balrog_output)
#    matchwithcosmos(ngmix_dfull, ngmix_dfull_match)
<<<<<<< HEAD
    stackwrite(y1a1_gold_tables, y1a1_gold)
#    matchwithcosmos(y1_d04, y1_d04_match_file)
#    chooseobjtype1(y1_d04)
#    chooseobjtype1(y1_d04_cosmos)
=======
#    stackwrite(y1a1_d04_tables, y1_d04)
    match(y1_d04_cosmos, y1_dfull_cosmos, 'd04', 'dfull', d04_dfull_match_file, d04_dfull)
#    matchwithcosmos(y1_d04, y1_d04_match_file)
>>>>>>> 0ba7e91eb5d38b3f66cef6d052ac3b389004dafc

def stackwrite(flist, output_file):
    tlist = [Table.read(table) for table in flist]
    new_table = vstack(tlist)
    del tlist
    new_table.write(output_file)

def match(file1, file2, name1, name2, match_file, outfile):
    tab1 = Table.read(file1)
    tab2 = Table.read(file2)

    h = esutil.htm.HTM(10)
    h.match(tab1['RA'], tab1['DEC'],
            tab2['RA'], tab2['DEC'],
            radius=1./3600,
            file=match_file)

    m = h.read(match_file)

    m1, m2, merr = np.array(zip(*m))

    print " {} matches".format(len(m1))
    print " No duplicates in matched list?:{}".format(len(m1)==len(set(m1)))

    m1_int = [int(g) for g in m1]
    m2_int = [int(c) for c in m2]

    new_tab = Table()
    for colname in tab1.colnames:
        new_tab.add_column(Column(data=tab1[colname][m1_int], name=colname+'_'+name1))
    for colname in tab2.colnames:
        new_tab.add_column(Column(data=tab2[colname][m2_int], name=colname+'_'+name2))

    new_tab.add_column(Column(name='match_err', data=merr))
    new_tab.write(outfile)

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



