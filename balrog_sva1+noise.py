import numpy as np
from astropy.table import Table

#sva1 gold with bad regions masks applied and cosmos matches removed (shouldn't matter)
sv = Table.read('/Users/Christina/DES/data/sva1_gold_auto_good_regions_no_cosmos.fits')
mc = (sv['MODEST_CLASS']==1)

#gather balrog sva1 tables
table_nums = []
for i in range(1, 12):
    if i==9:
        continue
    table_nums.append(str(i).zfill(2))

new_balrog = {}
for n in table_nums:
    b = Table.read('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'.format(n))
    for col in b.colnames:
        if n=='01':
            new_balrog[col] = b[col]
            new_balrog['TABLE_INDEX'] = ['01']
        else:
            new_balrog[col] = np.hstack([new_balrog[col], b[col]])
            new_balrog['TABLE_INDEX'].append(n)

#combined balrog table
balrog = Table(new_balrog)

#write combined table
balrog.write('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_SIM_TRUTH_zp_corr_fluxes.fits')

#reasonable bin cutoffs and sizes
bins = np.linspace(0, 500, 100)

for b in 'GRIZ':
    d = {}
    d['flux'] = sv['FLUX_AUTO_'+b][mc]
    d['err'] = sv['FLUXERR_AUTO_'+b][mc]

    d['flux0'] = balrog['FLUX_0_'+b]
    d['fluxNL'] = balrog['FLUX_NOISELESS_'+b]
    
    #bin the fluxes and errors
    sv_digitized = np.digitize(d['flux'], bins)
    b0_digitized = np.digitize(d['flux0'], bins)
    bNL_digitized = np.digitize(d['fluxNL'], bins)
    
    for ibin in len(1, bins):
        err_median = np.median(d['err'][sv_digitized == ibin])

        flux0_bin = d['flux0'][b0_digitized == ibin]
        fluxNL_bin = d['fluxNL'][bNL_digitized == ibin]

        errs_0 = np.random.normal(0., err_median, len(flux0_bin))
        errs_NL = np.random.normal(0., err_median, len(fluxNL_bin))

        flux0_bin = flux0_bin + errs_0
        fluxNL_bin = fluxNL_bin + errs_NL

        d['flux0'][b0_digitized == ibin] = flux0_bin
        d['fluxNL'][bNL_digitized == ibin] = fluxNL_bin

    balrog['FLUX_0+NOISE_'+b] = d['flux0']
    balrog['FLUX_NOISELESS+NOISE_'+b] = d['fluxNL']
    
balrog.write('/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_SIM_TRUTH_zp_corr_noised_fluxes.fits')
