import numpy as np
from astropy.table import Table, Column

sva1_gold = Table.read('/home/ckrawiec/DES/data/sva1_gold_detmodel_gals.fits')
cosmos = Table.read('/home/ckrawiec/COSMOS/data/COSMOS2015_Laigle+_v1.1.fits')

sv_no_cosmos = np.loadtxt('/home/ckrawiec/DES/magnification/lbgselect/gold_no_cosmos_indices.txt', dtype=int)
sv_cosmos = np.loadtxt('/home/ckrawiec/DES/magnification/lbgselect/gold_cosmos15_indices.txt', dtype=int)
cosmos_ind = np.loadtxt('/home/ckrawiec/DES/magnification/lbgselect/cosmos15_indices.txt', dtype=int)

new_sv = sva1_gold[sv_no_cosmos]
new_cosmos = sva1_gold[sv_cosmos] 
new_cosmos.add_column(Column(name='photoz', data=cosmos[cosmos_ind]))

new_sv.write('/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_no_cosmos.fits')
new_cosmos.write('/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_cosmos.fits')
