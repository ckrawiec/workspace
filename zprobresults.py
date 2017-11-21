import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from myutils import joincorrecttype
from scipy import stats
from scipy.misc import comb
from astropy.table import Table, join, Column

#uses z_column to check accuracy of results when target redshift is known
check_truth = False
z_column = ''
#'Z'
#'ZMINCHI2_dfull'
#'photoz_dfull'

#variables
num_files = 9
id_column_true = ''
#'BALROG_INDEX'
#'COADD_OBJECTS_ID_d04'
id_column_targ = 'COADD_OBJECTS_ID' 
#'BALROG_INDEX'
#'COADD_OBJECTS_ID_d04'
#'ID'
z_groups = [[0.1, 1.0],
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 9.9]]

id_col_type = float

truth_file = ''
#'/Users/Christina/DES/data/balrog/sva1/balrog_sva1_tab{}_TRUTH.fits'
#'/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits'
#'/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched_chunk{}-15.fits'
#'/Users/Christina/DES/data/y1a1_spec_gold_v2_weight_depth.fits'
data_file = ''
#'/Users/Christina/DES/data/y1a1_gold_flux_auto_griz_000001.fits'
#'/Users/Christina/DES/data/balrog/sva1/balrog_sva1_tab01_TRUTH.fits'
#'/Users/Christina/DES/data/balrog_sva1_TRUTH_zp_corr_fluxes_fixids.fits'
#'/Users/Christina/DES/data/y1a1_gold_d04_dfull_cosmos_matched.fits'
#'/Users/Christina/DES/data/sva1_gold_auto_good_regions_no_cosmos.fits'
#'/Users/Christina/DES/magnification/lbgselect/mocks/run7/zprob_mock{}.fits'

results_dir = '/home/ckrawiec/DES/magnification/lbgselect/zproboutput/'
#'/Users/Christina/DES/magnification/lbgselect/'
output_dir = '/home/ckrawiec/DES/magnification/lbgselect/zproboutput/'
name = 'y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_4bins_full_gauss_0000{}'
#'zprob_balrog_sva1_balrogz_z25_3bins_sigma_tree_noiseless_auto_both_zpcorr_griz_tab{}'
#'balrog_sva1_balrogz_noiseless_Efuncnoised_z3_3bins_full_gauss_tab01'
#'y1a1_gold_cosmosd04_cosmosdfull_zminchi2_auto_griz_full_z3_3bins'
#'y1a1_gold_cosmosd04_cosmosdfull_matched_zminchi2_auto_griz_full_2gauss_z3_4bins'
#'y1a1_gold_cosmosd04_cosmosdfull_matched_zminchi2_auto_griz_full_z3_4bins'
#'sva1_gold_cosmos_zminchi2_auto_griz_z3_3bins_full_gauss'
#'y1a1_gold_cosmosd04_cosmosdfull_matched_chunk{}-15_zminchi2_auto_griz_full_z0.5_3bins'
#'y1a1_spec_gold_v2_cosmosdfull_zminchi2_auto_griz_full_z3_3bins'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab{}'
#'y1a1_gold_cosmosdfull_zminchi2_auto_griz_z3_3bins_full_gauss_000001'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_f1s1_full_z3_3bins_Efunc_tab{}'
#'sva1_gold_cosmosdfull_photoz_auto_griz_z3_3bins_full_gauss_0-500000'
#'y1a1_gold_cosmosdfulld04noised_cosmosdfull_photoz_auto_griz_full_gauss_z3_2bins'
#'y1a1_gold_cosmosd04_cosmosdfull_matched_photoz_auto_griz_full_z3_3bins'
#'y1a1_gold_cosmosd04_cosmosdfull_zminchi2_auto_griz_full_gauss_z3_3bins'
#'y1a1_gold_cosmosd04_cosmosdfull_zminchi2_auto_griz_full_gauss_z3_3bins'
#'zprob_balrog_sva1_auto_cosmos_photoz_griz_zpcorr_full_z2.5_3bins_v2_tab{}'
#'balrog_sva1_balrogz_0griz_auto_zpcorr_f1s1_6culltree_z3_2bins_Efunc_tab{}'
#'run7_mock{}_newL_6culltree'

files = [''.join(os.path.splitext(f)[0]).split('/')[-1] for f in glob.glob(results_dir+name.format('*')+'*fits')]

def Pf(x, N, M):
    N, M = float(N), float(M)
    return comb(N, M) * x**M * (1.-x)**(N-M) * (N+1.)

def finderrors(xs, numlist, denlist, dx):
    errors = []
    for Mi, Ni in zip(numlist, denlist):
        frac = Mi/float(Ni)
        integ = np.array([np.sum(Pf(xs, Ni, Mi)[np.where((xs>=(frac-xi)) & (xs<=(frac+xi)))])*dx for xi in xs])
        error = xs[np.argmin(np.abs(integ-0.68))]
        errors.append(error)
    return errors

def autocolors(dtable, title, colormask): 
    gr = dtable['MAG_AUTO_G'] - dtable['MAG_AUTO_R']
    ri = dtable['MAG_AUTO_R'] - dtable['MAG_AUTO_I']

    grbox = ((gr < 4.) & (gr > -1.5))
    ribox = ((ri < 4.) & (ri > -1.5))

    plt.hist2d(ri[ribox & grbox & colormask], gr[ribox & grbox & colormask],
               bins=100)
    plt.xlabel('mag_auto_r - mag_auto_i')
    plt.ylabel('mag_auto_g - mag_auto_r')
    plt.colorbar()
    plt.savefig(title+'_autocolors.png')
    plt.close()

def main():
    print "Files found: "
    for fi in files:
        print "    ", fi
    
    #make file names
    nums = [str(i).zfill(2) for i in range(1,num_files+1)]
    names, trues = [],[]
    for f in files:
        if ('tab' in f) or ('chunk' in f) or ('_0000' in f):
            if f not in names:
                for num in nums:
                    if ('tab'+num in f) or ('chunk'+num in f) or ('_0000'+num in f):
                        names.append(f)
                        trues.append(truth_file.format(num))
                    
        else:
            names.append(f)
            trues.append(truth_file)
    
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
            results_dict['name'].append(names[i])
        else:
            print "File {} does not exist, moving to next one.".format(table_name)
            continue
    
        print "Working on ..."
        f.write("File group:\n")
        print "    ", table_name
        f.write("    {}\n".format(table_name))

        if check_truth==True:
            print "    ", trues[i]
            f.write("    {}\n".format(trues[i]))
        
            #join results and truth tables
            joined = joincorrecttype(trues[i], table_name,
                                     id_column_true, id_column_targ,
                                     id_col_type)
            f.write("Length of joined table: {}\n".format(len(joined)))
    
            #remove objects outside of redshift range
            joined = joined[(joined[z_column] >= (np.min(z_groups))) & (joined[z_column] < (np.max(z_groups)))]
            f.write("Length of joined table after truth redshift cuts: {}\n".format(len(joined)))
            f.write("\n"+names[i]+"\n")
    
        else:
            #joined = joincorrecttype(table_name, data_file,
            #                         id_column_targ, id_column_targ,
            #                         id_col_type)
            joined = Table.read(table_name)
            
            
        #loop through redshift groups and save P & N info
        for z_group in z_groups:
            P = joined['P'+str(z_group)]
    
            zp_mask = (P > 0.9)
            if 'MAG_AUTO_G' in joined.colnames:
                autocolors(joined, output_dir+names[i]+str(z_group), zp_mask)
            
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
    #set up colors for each z group (to loop through)
    colors = plt.cm.hsv(np.linspace(0, 0.8, len(z_groups)))

    ##z-prob histograms
    #for each file
    for i in range(len(results_dict['name'])):
        for z_group in z_groups:
            P = results_dict['P'+str(z_group)][i]
            plt.hist(P[~np.isnan(P)], bins=100,
                    histtype='step', log=True, label=str(z_group))
        plt.xlabel('P')
        plt.legend()
        plt.savefig(output_dir+results_dict['name'][i]+'_Phist.png')
        plt.close()
        
    #combined files
    for z_group in z_groups:
        P = np.hstack(results_dict['P'+str(z_group)])
        plt.hist(P[~np.isnan(P)], bins=50,
                 histtype='step', log=True, label=str(z_group))
        plt.xlabel('P')
    plt.legend()
    plt.savefig(output_dir+name.format('')+'_Phist.png')
    plt.close()
    
    ##saving probability information
    #probability bins
    bin_size = 0.05
    Pbins = np.arange(0., 1.0, bin_size)

    #bin midpoints
    xp = [ip+bin_size/2. for ip in Pbins]

    N_names = len(results_dict['name'])
 #   print N_names
 #   print results_dict['name']
    for z_group in z_groups:
        results_dict['single_Ns'+str(z_group)] = [[] for i in range(N_names)]
        results_dict['cut_above_Ns'+str(z_group)] = [[] for i in range(N_names)]
        results_dict['cut_below_Ns'+str(z_group)] = [[] for i in range(N_names)]
        results_dict['cut_above_sums'+str(z_group)] = [[] for i in range(N_names)]
        results_dict['cut_below_sums'+str(z_group)] = [[] for i in range(N_names)]
        results_dict['P_bin_sums'+str(z_group)] = [[] for i in range(N_names)]
        if check_truth==True:
            results_dict['correct'+str(z_group)] = [[] for i in range(N_names)]
#    print len(results_dict['single_Ns'+str(z_groups[0])])

    for i in range(len(results_dict['name'])):

        for z_group in z_groups:           
            P = results_dict['P'+str(z_group)][i]
            
            for ip in Pbins:
                #how many objects in each bin
                single_mask = (~np.isnan(P) & (P>ip) & (P<=(ip+bin_size)))
                #how many with P>Pbin
                cut_above_mask = (~np.isnan(P) & (P>=(ip)))
                #how many with P<Pbin
                cut_below_mask = (~np.isnan(P) & (P<(ip)))

                results_dict['single_Ns'+str(z_group)][i].append(len(P[single_mask]))
                results_dict['cut_above_Ns'+str(z_group)][i].append(len(P[cut_above_mask]))
                results_dict['cut_below_Ns'+str(z_group)][i].append(len(P[cut_below_mask]))
                results_dict['cut_above_sums'+str(z_group)][i].append(np.sum(P[cut_above_mask]))
                results_dict['cut_below_sums'+str(z_group)][i].append(np.sum(P[cut_below_mask]))
                results_dict['P_bin_sums'+str(z_group)][i].append(np.sum(P[single_mask]))
                if check_truth==True:
                    #number of objects in bin truly in z_group
                    results_dict['correct'+str(z_group)][i].append(float(len(P[results_dict['mask'+str(z_group)][i] & single_mask])))

    #print len(results_dict['single_Ns'+str(z_group)[0]])
    
    #plot them for each file
        
    for i in range(N_names):
        #bar hist
        ci=0
        for z_group in z_groups:
            plt.bar(xp, results_dict['single_Ns'+str(z_group)][i],
                    align='center', width=bin_size, color=colors[ci], alpha=0.2,
                    log=True, label=str(z_group))
#            plt.scatter(xp, results_dict['single_Ns'+str(z_group)][i],
#                        c=colors[ci])
            ci+=1
        plt.xlabel('P')
        plt.ylabel('N')
        plt.legend(loc='best')
        plt.savefig(output_dir+results_dict['name'][i]+'_Pbin_hist.png')
        plt.close()

        #cumulative totals
        ci=0
        for z_group in z_groups:
            plt.plot(xp, results_dict['cut_above_Ns'+str(z_group)][i],
                     'o-', markeredgecolor='none', c=colors[ci], label=str(z_group))
            plt.plot(xp, results_dict['cut_below_Ns'+str(z_group)][i],
                     'o--', markeredgecolor='none', c=colors[ci])
            ci+=1
        plt.legend(loc='best')
        plt.grid()
        plt.xlabel('P')
        plt.ylabel('cumulative number')
        plt.yscale('log')
        plt.savefig(output_dir+results_dict['name'][i]+'_Pbin_cumnumber.png')
        plt.close()

        #cumulative sums
        ci=0
        for z_group in z_groups:
            plt.plot(xp, results_dict['cut_above_sums'+str(z_group)][i],
                     'o-', markeredgecolor='none', c=colors[ci], label=str(z_group))
            plt.plot(xp, results_dict['cut_below_sums'+str(z_group)][i],
                     'o--', markeredgecolor='none', c=colors[ci])
            ci+=1
        plt.legend(loc='best')
        plt.grid()
        plt.xlabel('P')
        plt.ylabel('cumulative sums')
        plt.yscale('log')
        plt.savefig(output_dir+results_dict['name'][i]+'_Pbin_cumsums.png')
        plt.close()
    
    #combined files
    ci=0
    for z_group in z_groups:
        plt.bar(xp, np.sum(results_dict['single_Ns'+str(z_group)], axis=0),
                align='center', width=bin_size, color=colors[ci], alpha=0.2,
                log=True, label=str(z_group))
        plt.scatter(xp, np.sum(results_dict['single_Ns'+str(z_group)], axis=0),
                    edgecolor='none', c=colors[ci])
        ci+=1
    plt.xlabel('P')
    plt.ylabel('N')
    plt.legend(loc='best')
    plt.xlim(0., 1.)
    plt.savefig(output_dir+name.format('')+'_Pbin_hist.png')
    plt.close()

    ci=0
    for z_group in z_groups:
        plt.plot(xp, np.sum(results_dict['cut_above_Ns'+str(z_group)], axis=0),
                'o-', markeredgecolor='none', c=colors[ci], label=str(z_group))
        plt.plot(xp, np.sum(results_dict['cut_below_Ns'+str(z_group)], axis=0),
                'o--', markeredgecolor='none', c=colors[ci])
        ci+=1
    plt.xlabel('P')
    plt.ylabel('cumulative number')
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(output_dir+name.format('')+'_Pbin_cumnumber.png')
    plt.close()

    ci=0
    for z_group in z_groups:
        plt.plot(xp, np.sum(results_dict['cut_above_sums'+str(z_group)], axis=0),
                'o-', markeredgecolor='none', c=colors[ci], label=str(z_group))
        plt.plot(xp, np.sum(results_dict['cut_below_sums'+str(z_group)], axis=0),
                'o--', markeredgecolor='none', c=colors[ci])
        ci+=1
    plt.xlabel('P')
    plt.ylabel('cumulative sums')
    plt.grid()
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig(output_dir+name.format('')+'_Pbin_cumsums.png')
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
        xf = np.arange(0., 1., 0.001)
        dx = xf[1]-xf[0]
        #each file
        for i in range(len(results_dict['name'])):
            ci=0
            plt.plot(xp, xp, 'k--', lw=2.)
            for iz in range(len(z_groups)):
                Ms = np.array(all_tt[i][iz])
                Ns = np.array(all_lenP[i][iz])
                errors = finderrors(xf, Ms, Ns, dx)
                plt.errorbar(xp, Ms/Ns, yerr=errors, ecolor=colors[ci])
                plt.plot(xp, Ms/Ns, 'o-', c=colors[ci], label=str(z_groups[iz]))
                ci+=1
            plt.ylim(0., 1.)
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
            Ms = np.sum(all_tt, axis=0)[iz]
            Ns = np.sum(all_lenP, axis=0)[iz]
            errors = finderrors(xf, Ms, Ns, dx)
            plt.errorbar(xp, Ms/Ns, yerr=errors, ecolor=colors[ci])
            plt.plot(xp, Ms/Ns, 'o-', c=colors[ci], label=str(z_groups[iz]))
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

if __name__=="__main__":
    main()
    
