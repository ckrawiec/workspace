import numpy as np
from astropy.table import Table, join

nums = [str(ri).zfill(2) for ri in range(1,12)]
nums.remove('09')

bands = 'GRIZ'

truth_form = '/home/ckrawiec/DES/data/balrog_sva1_tab{}_TRUTH.fits'
sim_form = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM.fits'

def main():
    zpcorrect()

def simtruth():
    out_form = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM_TRUTH.fits'

    for num in nums:
        print "working on table "+num+"..."
        truth = Table.read(truth_form.format(num))
        sim = Table.read(sim_form.format(num))
        
        tab = join(truth, sim)
        
        tab.write(out_form.format(num))

def zpcorrect():
    #correct fluxes so they are on same zeropoint scale (30.0)

    out_form = '/home/ckrawiec/DES/data/balrog_sva1_auto_tab{}_SIM_zp_corr_fluxes.fits'

    for num in nums:
        print "working on table "+num+"..."
        sim = Table.read(sim_form.format(num))
        newsim = sim

        for band in bands:
            zp = np.array(np.abs(sim['MAG_AUTO_'+band] + 2.5*np.log10(sim['FLUX_AUTO_'+band])))
            flux_factor = 10**(-0.4 * zp) / 10**(-0.4 * 30.)
            newsim['FLUX_AUTO_'+band] = sim['FLUX_AUTO_'+band] * flux_factor
            newsim['FLUXERR_AUTO_'+band] = sim['FLUXERR_AUTO_'+band] * flux_factor
        
        newsim.write(out_form.format(num))

if __name__=="__main__":
	main()
