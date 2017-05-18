import os
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy import stats
from astropy.table import Table, join, Column

#uses z_column to check accuracy of results when target redshift is known
check_truth = True
z_column = 'Z'
#'photoz_dfull'

#variables
num_files = 12
id_column_true = 'BALROG_INDEX'
#'COADD_OBJECTS_ID'
id_column_targ = 'BALROG_INDEX'
#'COADD_OBJECTS_ID'
#'ID'
z_groups = [[0.001, 0.8],
            [0.8, 2.5],
            [2.5, 9.9]]

id_col_type = float

truth_file = '/Users/Christina/DES/data/balrog/sva1/balrog_sva1_tab{}_TRUTH.fits'
#'/Users/Christina/DES/data/balrog_sva1_TRUTH_zp_corr_fluxes_fixids.fits'
#'/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits
#'/Users/Christina/DES/magnification/lbgselect/mocks/run7/zprob_mock{}.fits'

results_dir = '/Users/Christina/DES/magnification/lbgselect/'
output_dir = '/Users/Christina/DES/magnification/lbgselect/zproboutput/'
name = 'zprob_balrog_sva1_balrogz_z25_3bins_sigma_tree_noiseless_auto_both_zpcorr_griz_tab{}'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab{}'
#'y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_3bins_full_gauss_000001'
#'sva1_gold_cosmos_zminchi2_auto_griz_z3_3bins_full_gauss'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_f1s1_full_z3_3bins_Efunc_tab{}'
#'sva1_gold_cosmosdfull_photoz_auto_griz_z3_3bins_full_gauss_0-500000'
#'y1a1_gold_cosmosdfulld04noised_cosmosdfull_photoz_auto_griz_full_gauss_z3_2bins'
#'y1a1_gold_cosmosd04_cosmosdfull_matched_photoz_auto_griz_full_z3_3bins'
#'y1a1_gold_cosmosd04_cosmosdfull_zminchi2_auto_griz_full_gauss_z3_3bins'
#'y1a1_gold_cosmosd04_cosmosdfull_zminchi2_auto_griz_full_gauss_z3_3bins'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab{}'
#'balrog_sva1_balrogz_noiseless_Efuncnoised_z3_3bins_full_Efunc_tab01'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_f1s1_6culltree_z3_2bins_Efunc_tab{}'
#'run7_mock{}_newL_6culltree'

files = [''.join(os.path.splitext(f)[0]).split('/')[-1] for f in glob.glob(results_dir+name.format('*')+'*fits')]

print "Files found: "
for fi in files:
    print "    ", fi

#make file names
nums = [str(i).zfill(2) for i in range(1,num_files+1)]
names, trues = [],[]
for f in files:
    for num in nums:
        if 'tab' in f:
            if 'tab'+num in f:
                names.append(f)
                trues.append(truth_file.format(num))
        else:
            names.append(f)
            trues.append(truth_file)

#print names, trues

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

    print "Working on ..."
    print "    ", table_name
    if check_truth==True:
        print "    ", trues[i]
    
        #join results and truth tables
        true = Table.read(trues[i])

        #this block is to fix tables where id_columns are in different formats
        keep = []
        for iid in range(len(tab[id_column_targ])):
            try:
                id_col_type(tab[id_column_targ][iid])
                keep.append(iid)
            except ValueError:
                continue

        tab = tab[keep]
        ids = tab[id_column_targ]
        tab.remove_column(id_column_targ)
        tab.add_column(Column(data=ids, name=id_column_targ, dtype=id_col_type))

        #print true[id_column_true][0], tab[id_column_targ][0]
        if (id_column_targ != id_column_true):
            tab.rename_column(id_column_targ, id_column_true)
        joined = join(true, tab, keys=id_column_true)
        f.write("Length of joined table: {}\n".format(len(joined)))

        #remove objects outside of redshift range
        joined = joined[(joined[z_column] > (np.min(z_groups)-0.001)) & (joined[z_column] < (np.max(z_groups)+0.001))]
        f.write("Length of joined table after truth redshift cuts: {}\n".format(len(joined)))
        f.write("\n"+names[i]+"\n")

    else:
        joined = tab
        
    #loop through redshift groups and save P & N info
    for z_group in z_groups:

        P = joined['P'+str(z_group)]

        results_dict['P'+str(z_group)].append(P)

        nans+=len(P[np.isnan(P)])

        if check_truth==True:
            z_mask = (joined[z_column] > np.min(z_group)) & (joined[z_column] < np.max(z_group))
            results_dict['mask'+str(z_group)].append(z_mask)
            z = joined[z_mask]

            results_dict['N'+str(z_group)].append([len(z)])

            f.write("    # True z in {}: {}\n".format(z_group, len(z)))

        f.write("    sum(P{}) = {}\n".format(z_group, np.sum(P[~np.isnan(P)])))
           
f.write("Number of NaNs: {}\n".format(nans))
print "results written to "+stats_file

###plotting
##histograms
#for each file
for i in range(len(results_dict['name'])):
    for z_group in z_groups:
        P = results_dict['P'+str(z_group)][i]
        plt.hist(P[~np.isnan(P)], bins=100,
                histtype='step', label=str(z_group))
    plt.xlabel('P')
    plt.legend()
    plt.savefig(output_dir+results_dict['name'][i]+'_Phist.png')
    plt.close()
        
#combined files
for z_group in z_groups:
        P = np.hstack(results_dict['P'+str(z_group)])
        plt.hist(P[~np.isnan(P)], bins=100,
                histtype='step', label=str(z_group))
plt.xlabel('P')
plt.legend()
plt.savefig(output_dir+name.format('')+'_Phist.png')
plt.close()

#probability bins
bin_size = 0.05
Pbins = np.arange(0., 1.0, bin_size)

#bin midpoints
xp = [ip+bin_size/2. for ip in Pbins]

#loop through colors for each z group
colors = plt.cm.brg(np.linspace(0, 1, len(z_groups)))

#for each file
tot_singles, tot_cut = [], []
for i in range(len(results_dict['name'])):   
    ci=0
    file_singles, file_cut = [], []
    for z_group in z_groups:
        single_Ns, cut_Ns = [], []
        P = results_dict['P'+str(z_group)][i]
        for ip in Pbins:
            #how many objects in each bin
            single_mask = (~np.isnan(P) & (P>ip) & (P<=(ip+bin_size)))
            #how many with P>Pbin
            cut_mask = (~np.isnan(P) & (P>=(ip)))
            single_Ns.append(len(P[single_mask]))
            cut_Ns.append(len(P[cut_mask]))
        plt.bar(xp, single_Ns, align='center',
                alpha=0.4, width=bin_size, color=colors[ci], log=True)
        plt.plot(xp, cut_Ns, 'o-', c=colors[ci], label=str(z_group))
        ci+=1
        file_singles.append(single_Ns)
        file_cut.append(cut_Ns)
    tot_singles.append(file_singles)
    tot_cut.append(file_cut)
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('P')
    plt.ylabel('N')
    plt.savefig(output_dir+results_dict['name'][i]+'_Pbin_totals.png')
    plt.close()

#combined files
ci=0
for iz in range(len(z_groups)):
    plt.bar(xp, np.sum(tot_singles, axis=0)[iz], align='center',
            width=bin_size, color=colors[ci], alpha=0.4, log=True)
    plt.plot(xp, np.sum(tot_cut, axis=0)[iz], 'o-', c=colors[ci], label=str(z_groups[iz]))
    ci+=1
plt.xlabel('P')
plt.ylabel('N')
plt.grid()
plt.legend(loc='best')
plt.savefig(output_dir+name.format('')+'_Pbin_totals.png')
plt.close()

                   
def checktruth():
    #bin midpoints
    xp = [ip+bin_size/2. for ip in Pbins]
    ##P sums    
    
    all_tt = []
    all_lenP = []
    #for each file
    for i in range(len(results_dict['name'])):
        ci = 0
        tt_groups = []
        lenPbin_groups = []
    
        #plot with 2 subplots, top twice as high
        gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])   
        for z_group in z_groups:
            P = results_dict['P'+str(z_group)][i]
            mt, tt, lenPbin = [],[],[]
            for ip in Pbins:
                mask = (~np.isnan(P) & (P>ip) & (P<=(ip+bin_size)))
                #sum of probabilities in bin
                meas_total = np.sum(P[mask])
                #number of objects in bin truly in z_group
                true_total = float(len(P[results_dict['mask'+str(z_group)][i] & mask]))
        
                mt.append(meas_total)
                tt.append(true_total)
                lenPbin.append(float(len(P[mask])))
            if np.any(mt):
                ax1.semilogy(xp, mt, 'o--', c=colors[ci], label=str(z_group))
                ax1.semilogy(xp, tt, 'o-', c=colors[ci])
                ax2.plot(xp, [0.]*len(xp), 'k--')
                ax2.plot(xp, np.array(mt)-np.array(tt), 'o-', c=colors[ci])
            
                ci+=1
            tt_groups.append(tt)
            lenPbin_groups.append(lenPbin)
    
        all_tt.append(tt_groups)
        all_lenP.append(lenPbin_groups)

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
    f.write("total sum of probabilities: {}\n".format(np.sum(Psums)))
    f.write("total number of objects: {}\n".format(np.sum(Nsums)))

if check_truth==True:
    checktruth()

f.close()
