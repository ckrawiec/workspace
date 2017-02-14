import numpy as np
from astropy.table import Table

#sva1 gold with bad regions masks applied and cosmos matches removed (shouldn't matter)
sv_file = '/home/ckrawiec/DES/data/sva1_gold_auto_good_regions_no_cosmos.fits'

balrog_combined_output = '/home/ckrawiec/DES/data/balrog_sva1_TRUTH_zp_corr_fluxes.fits'
balrog_noise_output = '/home/ckrawiec/DES/data/balrog_sva1_TRUTH_zp_corr_noised_fluxes.fits'

#gather balrog sva1 tables
table_nums = []
for i in range(1, 12):
    if i==9:
        continue
    table_nums.append(str(i).zfill(2))

new_balrog = {}
for n in table_nums:
    b = Table.read('/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'.format(n))

    if n=='01':
        new_balrog['BALROG_INDEX'] = [n+'_'+str(bi) for bi in b['BALROG_INDEX']]
        new_balrog['TABLE_INDEX'] = ['01']*len(b)
        for band in 'GRIZ':
            new_balrog['Z'] = b['Z']
            new_balrog['FLUX_0_'+band] = b['FLUX_0_'+band]
            new_balrog['FLUX_NOISED_'+band] = b['FLUX_NOISED_'+band]
            new_balrog['FLUX_NOISELESS_'+band] = b['FLUX_NOISELESS_'+band]

    else:
        for k in new_balrog.keys():
            if k=='TABLE_INDEX':
                new_balrog['TABLE_INDEX'] = np.hstack([new_balrog['TABLE_INDEX'], [n]*len(b)])                
            elif k=='BALROG_INDEX':
                new_balrog['BALROG_INDEX'] = np.hstack([new_balrog['BALROG_INDEX'], 
                                                       [n+'_'+str(bi) for bi in b['BALROG_INDEX']]])
            else:
                new_balrog[k] = np.hstack([new_balrog[k], b[k]])

del b

print "combined tables"

#combined balrog table
balrog = Table.read(balrog_combined_output)#new_balrog)

print "made new table"

#write combined table
#balrog.write(balrog_combined_output)

print "written to file"

#read in sv table, use galaxy-type objects
sv = Table.read(sv_file)
mc = (sv['MODEST_CLASS']==1)

print "sv data read in"

#reasonable bin cutoffs and sizes
bins = np.linspace(0, 500, 100)

for b in 'GRIZ':
    d = {}
    d['flux'] = sv['FLUX_AUTO_'+b][mc]
    d['err'] = sv['FLUXERR_AUTO_'+b][mc]

    d['flux0'] = balrog['FLUX_0_'+b]
    d['fluxNL'] = balrog['FLUX_NOISELESS_'+b]
    
    d['0err'] = np.array([np.nan]*len(balrog))
    d['NLerr'] = np.array([np.nan]*len(balrog))

    #bin the fluxes and errors
    sv_digitized = np.digitize(d['flux'], bins)
    b0_digitized = np.digitize(d['flux0'], bins)
    bNL_digitized = np.digitize(d['fluxNL'], bins)
    
    for ibin in range(1, len(bins)):
        err_median = np.median(d['err'][sv_digitized == ibin])

        flux0_bin = d['flux0'][b0_digitized == ibin]
        fluxNL_bin = d['fluxNL'][bNL_digitized == ibin]

        errs_0 = np.random.normal(0., err_median, len(flux0_bin))
        errs_NL = np.random.normal(0., err_median, len(fluxNL_bin))

        flux0_bin = flux0_bin + errs_0
        fluxNL_bin = fluxNL_bin + errs_NL

        d['0err'][b0_digitized == ibin] = np.abs(errs_0)
        d['NLerr'][bNL_digitized == ibin] = np.abs(errs_NL)
        d['flux0'][b0_digitized == ibin] = flux0_bin
        d['fluxNL'][bNL_digitized == ibin] = fluxNL_bin

    balrog['FLUX_0+NOISE_'+b] = d['flux0']
    balrog['FLUX_NOISELESS+NOISE_'+b] = d['fluxNL']
    balrog['FLUXERR_0+NOISE_'+b] = d['0err']
    balrog['FLUXERR_NOISELESS+NOISE_'+b] = d['NLerr']

del sv
    
print "writing balrog+noise to file..."
balrog.write(balrog_noise_output)

print "done"
