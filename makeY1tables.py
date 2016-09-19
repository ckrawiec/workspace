import numpy as np
import glob
from astropy.table import Table, vstack

y1a1_gold_tables = glob.glob('/home/ckrawiec/DES/data/y1a1_gold_flux_detmodel_MC1*')

tlist = []
for table in y1a1_gold_tables:
    tlist.append(Table.read(table))

new_Y1 = vstack(tlist)

new_Y1.write('/home/ckrawiec/DES/data/y1a1_gold_flux_detmodel_MC1.fits')

#cosmos = Table.read('/home/ckrawiec/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits')

#new_sv = sva1_gold[sv_no_cosmos]
#new_cosmos = sva1_gold[sv_cosmos] 
#new_cosmos.add_column(Column(name='photoz', data=cosmos['photoz'][cosmos_ind]))

#new_sv.write('/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_no_cosmos.fits')
#new_cosmos.write('/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_cosmos.fits')
