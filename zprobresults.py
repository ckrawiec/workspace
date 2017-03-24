import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats
from astropy.table import Table, join, Column

#variables
output_dir = '/Users/Christina/DES/magnification/lbgselect/zproboutput/'
num_files = 6
id_column = 'BALROG_INDEX'

truth_file = '/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM_TRUTH_zp_corr_fluxes.fits'
#'/Users/Christina/DES/magnification/lbgselect/mocks/run1/zprob_mock{}.fits'

results_dir = '/Users/Christina/DES/magnification/lbgselect/zproboutput/'
name = 'balrog_sva1_balrogz_0griz_auto_zpcorr_5culltree_z3_2bins_tab{}_0-5000'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_5culltree_z3_2bins_Efunc_tab{}_0-10000'
#'balrog_sva1_flux0_SVnoised_griz_z3_2bins_4tree_tab{}_0-100000_culltree'
#'balrog_sva1_flux0_SVnoised_griz_z3_2bins_4tree_tab{}_0-10000_culltree_v2'
#'run1_mock{}_kmeans_tree_results_timetest_new'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_5tree_z3_3bins_tab{}_0-10000'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab{}'
#'run5_mock{}_4tree'
#'zprob_balrog_sva1_balrogz_z25_3bins_sigma_tree_noiseless_auto_both_zpcorr_griz_tab{}'
#'balrog_sva1_flux0_flat10noised_griz_z3_3bins_4tree_tab{}'
#'zprob_balrog_sva1_balrogz_auto_griz_zpcorr_0flat10noised_tree20k_z2.5_3bins_v2_tab{}'
#'zprob_balrog_sva1_balrogz_griz_zpcorr_tree_z2.5_3bins_v2_tab'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab'

z_column = 'Z'
z_groups = [[0.001, 3.0],
        #    [1.0, 3.0],
            [3.0, 9.9]]

#make file names
nums = [str(i).zfill(2) for i in range(1,num_files+1)]
names = [name.format(num) for num in nums]

trues = [truth_file.format(num) for num in nums]

#open text file
stats_file = output_dir+name.format('')+'_Pstats.txt'
f = open(stats_file, 'w')

#keep track of NaNs
nans = 0

#dict to save results from all tables
results_dict = {}
for z_group in z_groups:
    results_dict['P'+str(z_group)] = []
    results_dict['N'+str(z_group)] = []
    results_dict['mask'+str(z_group)] = []
    results_dict['name'] = []
    
#loop over results tables
for i in range(len(names)):
    table_name = results_dir+names[i]+'.fits'
    if os.path.exists(table_name):
        tab = Table.read(table_name)
        results_dict['name'].append(names[i])
    else:
        print "File {} does not exist, moving to next one.".format(table_name)
        continue

    #join results and truth tables
    true = Table.read(trues[i])

    keep = []
    for id in range(len(tab[id_column])):
        try:
            float(tab[id_column][id])
            keep.append(id)
        except ValueError:
            continue

    tab = tab[keep]
    ids = tab[id_column]
    tab.remove_column(id_column)
    tab.add_column(Column(data=ids, name=id_column, dtype='float'))

    joined = join(true, tab)
      
    f.write("\n"+names[i]+"\n")

    #loop through redshift groups and save P & N info
    for z_group in z_groups:

        P = joined['P'+str(z_group)]
        results_dict['P'+str(z_group)].append(P)

        nans+=len(P[np.isnan(P)])
        
        z_mask = (joined['Z'] > np.min(z_group)) & (joined['Z'] < np.max(z_group))
        results_dict['mask'+str(z_group)].append(z_mask)
        z = joined[z_mask]
        results_dict['N'+str(z_group)].append([len(z)])

        f.write("    # True z in {}: {}\n".format(z_group, len(z)))
        f.write("    sum(P{}) = {}\n".format(z_group, np.sum(P[~np.isnan(P)])))
           
print "results written to "+stats_file
print "N NaNs: ", nans
f.close()

###plotting
##histograms
#for each file
for i in range(len(results_dict['name'])):
    for z_group in z_groups:
        P = results_dict['P'+str(z_group)][i]
        plt.hist(P[~np.isnan(P)], histtype='step', label=str(z_group))
    plt.xlabel('P')
    plt.legend()
    plt.savefig(output_dir+results_dict['name'][i]+'_Phist.png')
    plt.close()
    
#combined files
P = np.hstack(results_dict['P'+str(z_group)])
for z_group in z_groups:
    plt.hist(P[~np.isnan(P)], histtype='step', label=str(z_group))
plt.xlabel('P')
plt.legend()
plt.savefig(output_dir+name.format('')+'_Phist.png')
plt.close()


##P sums    
#loop through colors for each z group
colors = plt.cm.brg(np.linspace(0, 1, len(z_groups)))
#probability bins
bin_size = 0.05
Pbins = np.arange(0., 1.0, bin_size)

all_tt = []
all_lenP = []
#for each file
for i in range(len(results_dict['name'])):
    ci = 0
    tt_groups = []
    lenP_groups = []
    #plot with 2 subplots, top twice as high
    gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])   
    for z_group in z_groups:
        P = results_dict['P'+str(z_group)][i]
        xp, mt, tt, lenP = [],[],[],[]
        for ip in Pbins:
            xp.append(ip+bin_size/2.)
            mask = (~np.isnan(P) & (P>ip) & (P<=(ip+bin_size)))
            meas_total = np.sum(P[mask])
            true_total = float(len(P[results_dict['mask'+str(z_group)][i] & mask]))

            mt.append(meas_total)
            tt.append(true_total)
            lenP.append(float(len(P[mask])))

        ax1.semilogy(xp, mt, 'o--', c=colors[ci], label=str(z_group))
        ax1.semilogy(xp, tt, 'o-', c=colors[ci])
        ax2.plot(xp, [0.]*len(xp), 'k--')
        ax2.plot(xp, np.array(mt)-np.array(tt), 'o-', c=colors[ci])
        
        ci+=1
        tt_groups.append(tt)
        lenP_groups.append(lenP)
    all_tt.append(tt_groups)
    all_lenP.append(lenP_groups)
            
    ax1.legend(loc='upper left')
    ax1.set_ylabel('Sum(P) (--) or # True (-)')
    ax1.grid()
    ax2.grid()
    ax2.set_xlabel('P')
    plt.subplots_adjust(hspace=0.)
    plt.ylabel('Sum(P) - True')
    plt.savefig(output_dir+results_dict['name'][i]+'_Prange_totals.png')
    plt.close()

#all files
#plot with 2 subplots, top twice as high
gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])   
ci=0
for z_group in z_groups:
    P = np.hstack(results_dict['P'+str(z_group)])

    xp, mt, tt = [],[],[]
    for ip in Pbins:
        xp.append(ip+bin_size/2.)
        mask = (~np.isnan(P) & (P>ip) & (P<=(ip+bin_size)))
        meas_total = np.sum(P[mask])
        z_mask = np.hstack(results_dict['mask'+str(z_group)])
        true_total = len(P[z_mask & mask])

        mt.append(meas_total)
        tt.append(true_total)

    ax1.semilogy(xp, mt, 'o--', c=colors[ci], label=str(z_group))
    ax1.semilogy(xp, tt, 'o-', c=colors[ci])
    ax2.plot(xp, [0.]*len(xp), 'k--')
    ax2.plot(xp, np.array(mt)-np.array(tt), 'o-', c=colors[ci])
        
    ci+=1
        
ax1.legend(loc='upper left')
ax1.set_ylabel('Sum(P) (--) or # True (-)')
ax1.grid()
ax2.grid()
ax2.set_xlabel('P')
plt.subplots_adjust(hspace=0.)
plt.ylabel('Sum(P) - True')
plt.savefig(output_dir+name.format('')+'_Prange_totals.png')
plt.close()

##true ratio in P bins
#each file
for i in range(len(results_dict['name'])):
    ci=0
    plt.plot(xp, xp, 'k--', lw=2.)
    for iz in range(len(z_groups)):
        plt.plot(xp, np.array(all_tt[i][iz])/np.array(all_lenP[i][iz]), 'o-', c=colors[ci], label=str(z_groups[iz]))
        ci+=1
    plt.xlabel('P')
    plt.ylabel('Fraction of objects truly in redshift range')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(output_dir+results_dict['name'][i]+'_Prange_trueratio.png')
    plt.close()
    
#combined files
ci=0
plt.plot(xp, xp, 'k--', lw=2.)
for iz in range(len(z_groups)):
    plt.plot(xp, np.sum(all_tt, axis=0)[iz]/np.sum(all_lenP, axis=0)[iz], 'o-', c=colors[ci], label=str(z_groups[iz]))
    ci+=1
plt.xlabel('P')
plt.ylabel('Fraction of objects truly in redshift range')
plt.grid()
plt.legend(loc='best')
plt.savefig(output_dir+name.format('')+'_Prange_trueratio.png')
plt.close()

#plot fractional difference between sum(P) & N
x = range(1, len(results_dict['name'])+1)
plt.plot(x, [0]*len(x), c='k')

Psums = []
Nsums = []
for i in range(len(results_dict['name'])):
    for z_group in z_groups:
        P = results_dict['P'+str(z_group)][i]
        Psum = np.sum(P[~np.isnan(P)])
        Psums.append(Psum)
        Nsum = np.sum(results_dict['N'+str(z_group)][i])
        Nsums.append(Nsum)
        plt.plot(x[i], (Psum-Nsum)/Nsum, 'o', label=str(z_group))

plt.ylabel('fractional difference')
plt.xlabel('file')
plt.title('Sum(P) vs. True N')
plt.legend()
plt.grid()
plt.savefig(output_dir+name.format('')+'_Psums.png')
plt.close()
        
#print totals
print "total Psum: ", np.sum(Psums)
print "N total: ", np.sum(Nsums)
