import numpy as np
from astropy.table import Table, join

nums = [str(ri).zfill(2) for ri in range(1,12)]
nums.remove('09')

truth_form = '/home/ckrawiec/DES/data/balrog_sva1_tab{}_TRUTH.fits'
sim_form = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM.fits'
out_form = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM_TRUTH.fits'

for num in nums:
    print "working on table "+num+"..."
    truth = Table.read(truth_form.format(num))
    sim = Table.read(sim_form.format(num))

    tab = join(truth, sim)
    
    tab.write(out_form.format(num))

