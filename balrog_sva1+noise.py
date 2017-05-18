import numpy as np
import sys
from astropy.table import Table

error_type = 'y1d04'  #'flat'
flat_error_sigma = 10.

#sva1 gold with bad regions masks applied and cosmos matches removed (shouldn't matter)
sv_file = '/Users/Christina/DES/data/sva1_gold_auto_good_regions_no_cosmos.fits'

#y1a1 gold all d04 objects that match to dfull-cosmos15 matched catalog
y1d04_file = '/Users/Christina/DES/data/y1a1_gold_d04_cosmos_matched.fits'

#balrog_combined_output = '/Users/Christina/DES/data/balrog_sva1_TRUTH_zp_corr_fluxes.fits'
balrog_noise_output = '/Users/Christina/DES/data/balrog_sva1_tab{}_TRUTH_zp_corr_SVnoised_fluxes.fits'

Efunc_params = [  0.3434363,
                -0.20932157,
                -0.13738978,
                 0.80368388,
                 4.48815913,
                 2.64501352]

def Gfunc(x, A, C, V):
    return A * np.exp(-0.5 * (x-C)*(x-C) / V)
    
def Efunc(x, A, B, C1, C2, V1, V2):
    return Gfunc(x, A, C1, V1) + Gfunc(x, B, C2, V2)

def eflat(flux, background):
    sigma = background
    return np.random.normal(0., sigma, len(flux)), [sigma]*len(flux)

values = np.arange(-100, 100, 0.001)
Eargs = [values]+Efunc_params
ps = Efunc(*Eargs)
ps /= sum(ps)

def enew(flux, background):
    sigmas = np.sqrt(flux + background**2.)
    return np.prod(zip(np.random.choice(values, len(flux), p=ps), sigmas), axis=1), sigmas


#gather balrog sva1 tables
table_nums = []
for i in range(1, 12):
    if i==9:
        continue
    table_nums.append(str(i).zfill(2))

if error_type == 'sv':
    #read in sv table, use galaxy-type objects
    sv = Table.read(sv_file)
    mc = (sv['MODEST_CLASS']==1)
    print "sv data read in"

    #reasonable bin cutoffs and sizes
    bins = np.linspace(0, 500, 100)

balrog = {}
for n in table_nums:
    btab = Table.read('/Users/Christina/DES/data/balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'.format(n))

    balrog['BALROG_INDEX'] = btab['BALROG_INDEX']
    balrog['TABLE_INDEX'] = [n]*len(btab)
    balrog['Z'] = btab['Z']

    for b in 'GRIZ':
        balrog['FLUX_0_'+b] = btab['FLUX_0_'+b]
        balrog['FLUX_NOISED_'+b] = btab['FLUX_NOISED_'+b]
        balrog['FLUX_NOISELESS_'+b] = btab['FLUX_NOISELESS_'+b]
                
        if error_type == 'sv':
            d = {}
            d['flux'] = sv['FLUX_AUTO_'+b][mc]
            d['err'] = sv['FLUXERR_AUTO_'+b][mc]
            
            d['flux0'] = balrog['FLUX_0_'+b]
            d['fluxNL'] = balrog['FLUX_NOISELESS_'+b]
            
            d['0err'] = np.array([np.nan]*len(d['flux0']))
            d['NLerr'] = np.array([np.nan]*len(d['fluxNL']))
            
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

               d['0err'][b0_digitized == ibin] = [err_median]*len(flux0_bin)
               d['NLerr'][bNL_digitized == ibin] = [err_median]*len(fluxNL_bin)
       ###THIS IS CREATING PROBLEM WHERE ORIGINAL FLUX0 COLUMN ALSO CHANGES TO NEW NOISED VALUE (MAKE EXPLICIT COPY?)        d['flux0'][b0_digitized == ibin] = flux0_bin
               d['fluxNL'][bNL_digitized == ibin] = fluxNL_bin

            balrog['FLUX_0+NOISE_'+b] = d['flux0']
            balrog['FLUX_NOISELESS+NOISE_'+b] = d['fluxNL']
            balrog['FLUXERR_0+NOISE_'+b] = d['0err']
            balrog['FLUXERR_NOISELESS+NOISE_'+b] = d['NLerr']

        elif error_type == 'flat':
            errors = np.random.normal(0., flat_error_sigma, len(balrog['FLUX_0_G']))
            balrog['FLUX_0+NOISE_'+b] = balrog['FLUX_0_'+b] + errors
            balrog['FLUX_NOISELESS+NOISE_'+b] = balrog['FLUX_NOISELESS_'+b] + errors
            balrog['FLUXERR_0+NOISE_'+b] = [flat_error_sigma]*len(errors)
            balrog['FLUXERR_NOISELESS+NOISE_'+b] = [flat_error_sigma]*len(errors)

        else:
            sys.exit("Choose 'flat' or 'sv' errors. Exiting")

    print "writing tab{} to file...".format(n)
    balrog_tab = Table(balrog)
    balrog_tab.write(balrog_noise_output.format(n))

print "done"
