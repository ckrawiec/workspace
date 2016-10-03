import numpy as np
import glob
import esutil
from astropy.table import Table, Column, vstack

home_dir = '/home/ckrawiec/'

def stackwrite(flist, output_file):
    tlist = [Table.read(table) for table in flist]
    new_table = vstack(tlist)
    del tlist
    new_table.write(output_file)

#*********BALROG Y1A1 TRUTH/SIM INNER JOIN, DETMODEL FLUXES********    
#balrog_tables = glob.glob('/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel_*')
#stackwrite(balrog_tables, '/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel.fits')

#*********Y1A1 GOLD DFULL********
#y1a1_gold_tables = glob.glob('/home/ckrawiec/DES/data/y1a1_gold_flux_detmodel_MC1_*')
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


#*********Y1A1 GOLD DFULL/COSMOS MATCHED BY POSITION********
cosmos = Table.read(home_dir+'COSMOS/data/COSMOS2015_Laigle+_v1.1.fits')

y1 = Table.read(home_dir+'DES/data/y1a1_gold_dfull.fits')

h = esutil.htm.HTM(10)
h.match(y1['RA'], y1['DEC'],
        cosmos['ALPHA_J2000'], cosmos['DELTA_J2000'],
        radius=1./3600,
        file=home_dir+'DES/data/match_y1a1_gold_dfull_cosmos_1arcsec')
m = h.read(home_dir+'DES/data/match_y1a1_gold_dfull_cosmos_1arcsec')

gold_m, cosmos_m, merr = np.array(zip(*m))

print "{} COSMOS/Y1 matches".format(len(cosmos_m))

print "matched"

no_cosmos = list(set(range(0,len(y1)))-set(gold_m))

print "No duplicates in matched list?:{}".format(len(gold_m)==len(set(gold_m)))

new_y1 = y1[no_cosmos]

print "made new Y1 table"

new_y1.write(home_dir+'DES/data/y1a1_gold_dfull_no_cosmos.fits')
del new_y1

print "wrote new Y1 table"

gold_m_int = [int(g) for g in gold_m]
cosmos_m_int = [int(c) for c in cosmos_m]

new_cosmos = y1[gold_m_int] 
del y1

print "deleted old Y1 table"

new_cosmos.add_column(Column(name='photoz', data=cosmos['PHOTOZ'][cosmos_m_int]))
new_cosmos.add_column(Column(name='NUMBER', data=cosmos['NUMBER'][cosmos_m_int]))
new_cosmos.add_column(Column(name='ALPHA_J2000', data=cosmos['ALPHA_J2000'][cosmos_m_int]))
new_cosmos.add_column(Column(name='DELTA_J2000', data=cosmos['ALPHA_J2000'][cosmos_m_int]))
new_cosmos.add_column(Column(name='match_err', data=merr))

print "made new cosmos table"

new_cosmos.write(home_dir+'DES/data/y1a1_gold_dfull_cosmos.fits')
