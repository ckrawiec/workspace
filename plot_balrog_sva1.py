import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import ckdtree
from astropy.table import Table, join
from collections import defaultdict

base_name = '/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_auto_cosmos_photoz_griz_full_z2.5_3bins_tab01_v2'
#'zprob_balrog_sva1_z25_3bins_photoz_sigma_tree_auto_zpcorr_griz'
#zprob_balrog_sva1_balrogz_z25_3bins_sigma_tree_noiseless_auto_both_zpcorr_griz'
#zprob_balrog_sva1_balrogz_z25_3bins_sigma_tree_auto_zpcorr_griz'

cosmos_file = '/home/ckrawiec/DES/data/sva1_gold_auto_good_regions_cosmos.fits'
#'/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_good_regions_cosmos.fits'

z_col = 'zminchi2'
SE_col = 'AUTO'
filters = ['G','R','I','Z']

zprob_file = '{}_tab{}.fits'
sim_file = '/home/ckrawiec/DES/data/balrog_sva1_'+SE_col.lower()+'_tab{}_SIM.fits'
truth_file = '/home/ckrawiec/DES/data/balrog_sva1_tab{}_TRUTH.fits'

num_list = [str(i).zfill(2) for i in range(1,12)]
num_list.remove('09')

#z_groups = [[0.001, 0.5],
#            [0.5, 1.0],
#            [1.0, 2.0],
#            [2.0, 3.0],
#            [3.0, 4.0],
#            [4.0, 5.0],
#            [5.0, 9.9]]
z_groups = [[0.001, 0.8],
            [0.8, 2.5],
            [2.5, 9.9]]

Pdict = {}
for z_group in z_groups:
    Pdict[str(z_group)] = []

truthz, ra, dec = [],[],[]
g,r,i,z = [],[],[],[]
gmag,rmag,imag,zmag = [],[],[],[]
gerr,rerr,ierr,zerr = [],[],[],[]
idmask, flags, objtype = [],[],[]

def main():
    #plotStats()
    printStats(0.7)
    #radec()
    #Phightruthz()
    #magi(cut=0.7)
    #truthzhist(0.7)
    #checkImposters(cut=0.8)

def list_uniques(seq):
    tally = defaultdict(list)
    for it, item in enumerate(seq):
        tally[item].append(it)
    return np.hstack((locs for key,locs in tally.items() if len(locs)==1))

for num in num_list:
    zprob = Table.read(zprob_file.format(base_name,num))
    sim = Table.read(sim_file.format(num))
    truth = Table.read(truth_file.format(num))

    tab1 = join(zprob, sim)
    tab = join(tab1, truth)

    for z_group in z_groups:

        Pdict[str(z_group)].append(tab['P'+str(z_group)])

    #remove duplicate ids
    falsemask = np.array([False]*len(tab))
    falsemask[list_uniques(tab['BALROG_INDEX'])]=True
    idmask.append(falsemask)

    flags.append(np.array(zip(*[tab['FLAGS_'+band] for band in filters])))
    objtype.append(tab['OBJTYPE'])
    truthz.append(tab['Z'])

    ra.append(tab['RA'])
    dec.append(tab['DEC'])
    
    g.append(tab['FLUX_'+SE_col+'_G'])
    r.append(tab['FLUX_'+SE_col+'_R'])
    i.append(tab['FLUX_'+SE_col+'_I'])
    z.append(tab['FLUX_'+SE_col+'_Z'])
    
    imag.append(tab['MAG_'+SE_col+'_I'])

    gerr.append(tab['FLUXERR_'+SE_col+'_G'])
    rerr.append(tab['FLUXERR_'+SE_col+'_R'])
    ierr.append(tab['FLUXERR_'+SE_col+'_I'])
    zerr.append(tab['FLUXERR_'+SE_col+'_Z'])

flags = np.vstack(flags)
objtype = np.hstack(objtype)

truthz = np.hstack(truthz)
ra = np.hstack(ra)
dec = np.hstack(dec)

g = np.hstack(g)
r = np.hstack(r)
i = np.hstack(i)
z = np.hstack(z)

imag = np.hstack(imag)

gerr = np.hstack(gerr)
rerr = np.hstack(rerr)
ierr = np.hstack(ierr)
zerr = np.hstack(zerr)

for k in Pdict.keys():
    Pdict[k] = np.hstack(Pdict[k])

#MASKS
idmask = np.hstack(idmask)
#remove objects with any flags
flagmask = np.array([np.all(iflag) for iflag in (flags==0)])
galmask = (objtype==1)

print base_name.split('/')[-1]
print "# Total Balrog objects: "+str(len(g))
print "# Balrog objects with unique ids: "+str(len(g[idmask]))
print "# Balrog objects with no flags: "+str(len(g[flagmask]))
print "# Balrog objects with OBJTYPE=1 (galaxies): "+str(len(g[galmask]))
print "# Balrog objects with unique ids, no flags, and OBJTYPE=1: \n"+str(len(g[idmask & flagmask & galmask]))

def plotStats():
    for z_group in z_groups:
        truthzCutHist(Pdict[str(z_group)], z_group)

def truthzCutHist(Pgrp, zgrp):
    zmask = (truthz>np.min(zgrp)) & (truthz<np.max(zgrp))
    for cut in np.arange(0,1.0,0.2):
        plt.hist(truthz[zmask & (Pgrp>cut)], histtype='step', label='P({})>{}'.format(zgrp, cut))
    plt.xlabel('truth z')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('{}_truthz_{}cut_hist.png'.format(base_name, zgrp))
    plt.close()

    for cut in np.arange(0,1.0,0.2):
        plt.hist(truthz[Pgrp>cut], histtype='step', label='P({})>{}'.format(zgrp, cut))
    plt.xlabel('truth z')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('{}_truthz_{}cut_hist_allz.png'.format(base_name, zgrp))
    plt.close()

def magi(cut):
    Phigh = Pdict[str(z_groups[-1])]
    plt.hist(imag[Phigh>cut], histtype='step', bins=1000, label='all')
    plt.yscale('log')
    plt.xlim(17,26)
    plt.legend(loc='upper left')
    plt.xlabel('mag_i')
    plt.title('P{} > {}'.format(z_groups[-1], cut))
    plt.savefig('{}_hist_mag_i_Pcut_high'.format(base_name))
    plt.close()

def printStats(cut):
    fig, ax = plt.subplots()
    fig.canvas.draw()

    for z_group in z_groups:
        truthzmask = (truthz>=(np.min(z_group))) & (truthz<(np.max(z_group)))
        P = Pdict[str(z_group)]
        
        purity, completeness, labels = [],[],[]

        def withmask(mask, title):
	    NP = len(truthz[mask & (P>cut)])        
            Ntrue = len(truthz[mask & truthzmask])
            Ncorrect = len(truthz[mask & truthzmask  & (P>cut)])
            print title 
            print "Number of galaxies with P{} > {}: {}".format(z_group, 
                                                                cut,
                                                                NP)
            print "Number of galaxies with truth z in  {}: {}".format(z_group, Ntrue)
            print "Number of galaxies with truth z in {} & P{} > {}: {} \n\
            ({}% complete, {}% pure)\n".format(z_group, 
                                               z_group, 
                                               cut, 
                                               Ncorrect, 
                                               Ncorrect/float(Ntrue) * 100.,
                                               Ncorrect/float(NP) * 100.)
            
            labels.append(title)
            purity.append(Ncorrect/float(NP) * 100.)
            completeness.append(Ncorrect/float(Ntrue) * 100.)

        withmask(mask=(np.array([True]*len(truthz))), title='#ALL objects')
        withmask(mask=(idmask), title='#Objects with NO duplicate ids')
        withmask(mask=(flagmask), title='#Objects with NO flags')
        withmask(mask=(galmask), title='#Objects with OBJTYPE=1')
        withmask(mask=(idmask & flagmask), title='#Objects with NO duplicate ids and NO flags')
        withmask(mask=(flagmask & galmask), title='#Objects with NO flags and OBJTYPE=1')
        withmask(mask=(idmask & flagmask & galmask), title='#Objects with NO duplicate ids, NO flags, and OBJTYPE=1')

        
        x = range(7)
        plt.plot(x, purity, 'o-', label=str(z_group))
        plt.plot(x, completeness, 'o--', label=str(z_group))

    pltlabels = [item.get_text() for item in ax.get_xticklabels()]
    for i in range(len(pltlabels)):
        pltlabels[i] = labels[i]

    ax.set_xticklabels(pltlabels)
    plt.legend(loc='best')
    plt.ylabel('purity(-)/completeness(--)')
    plt.savefig('{}_stats'.format(base_name))
    plt.close()

def radec():    
    plt.scatter(np.hstack(ra), np.hstack(dec), edgecolor='none')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.savefig('/home/ckrawiec/DES/data/balrog_sva1_SIM_ra_dec')
    plt.close()
    
def Phightruthz():
    Phigh = Pdict[str(z_groups[-1])]
    Pmid = Pdict[str(z_groups[1])]
    Plow = Pdict[str(z_groups[0])]
    plt.scatter(truthz, Phigh,
                s=4., edgecolor='none', c='b', label=str(z_groups[-1]))
    plt.scatter(truthz, Plow,
                s=4., edgecolor='none', c='g', label=str(z_groups[0]))
    plt.scatter(truthz, Pmid,
                s=4., edgecolor='none', c='c', label=str(z_groups[1]))
    plt.legend(loc='best')
    plt.xlabel('truth z')
    plt.ylabel('P')
    plt.savefig('{}_truthz_P_scatter'.format(base_name))
    plt.close()

    notnan = (~np.isnan(Plow)) & (~np.isnan(Phigh))

    plt.hist2d(truthz[notnan], Phigh[notnan], bins=100.)
    plt.xlabel('truth z')
    plt.ylabel('P[2.5, 9.9]')
    plt.savefig('{}_truthz_Phigh_hist2d'.format(base_name))
    plt.close()

def truthzhist(cut):
    for z_group in z_groups:
        plt.hist(truthz[Pdict[str(z_group)]>cut], histtype='step',label='P{} > {}'.format(z_group, cut), bins=100)
    plt.yscale('log')
    plt.xlabel('true z')
    plt.legend(loc='best')
    plt.savefig('{}_z_hist'.format(base_name))
    plt.close()

def checkImposters(cut):
    cosmos = Table.read(cosmos_file)

    Phigh = Pdict[str(z_groups[-1])]
    Plow = Pdict[str(z_groups[0])]

    wrong = (Phigh>cut) & (truthz<np.min(zsrc)) 
    right = (Phigh>cut) & (truthz>np.min(zsrc)) & (truthz<np.max(zsrc))
    nfilt = range(4)
    ncheck = 10
    
    z_mask = (cosmos[z_col] >= np.min(zsrc)) & (cosmos[z_col] < np.max(zsrc))

    data_zip = np.array( zip( *[g,r,i,z] ) )
    err_zip = np.array( zip( *[gerr, rerr, ierr, zerr] ) )

    template_zip = np.array( zip( *[cosmos['FLUX_'+SE_col+'_'+f][z_mask] for f in filters] ) )
    sigma = np.median(err_zip.T, axis=1)

    cosmos_gmag = cosmos['MAG_'+SE_col+'_G']
    cosmos_rmag = cosmos['MAG_'+SE_col+'_R']
    cosmos_imag = cosmos['MAG_'+SE_col+'_I']
    cosmos_zmag = cosmos['MAG_'+SE_col+'_Z']

    print "working on tree..."
    truetree = ckdtree.cKDTree(template_zip/sigma)
    print "done building tree"

    for icheck in range(ncheck):
        print "checking imposter #{}".format(icheck)
        dnear, inear = truetree.query(data_zip[wrong][icheck]/sigma, k=15)
        truearr = truetree.data[inear] * sigma
        
        for truearri in truearr:
            plt.plot(nfilt, truearri, 'o-', c='b')
        plt.plot(nfilt, data_zip[wrong][icheck], 'o-', c='g')
        plt.xlabel('filter (g,r,i,z)')
        plt.ylabel('flux_detmodel')
        plt.title('truth z = {}, Phigh = {}'.format(truthz[wrong][icheck], Phigh[wrong][icheck]))
        plt.savefig('{}_fluxes_imposter_{}'.format(base_name, str(icheck).zfill(3)))
        plt.close()
        
        for template in template_zip:
            plt.plot(nfilt, template, 'o-', c='r')
        for truearri in truearr:
            plt.plot(nfilt, truearri, 'o-', c='b')
        plt.plot(nfilt, data_zip[wrong][icheck], 'o-', c='g')
        plt.xlabel('filter (g,r,i,z)')
        plt.ylabel('flux_detmodel')
        plt.title('truth z = {}, Phigh = {}'.format(truthz[wrong][icheck], Phigh[wrong][icheck]))
        plt.savefig('{}_all_fluxes_imposter_{}'.format(base_name, str(icheck).zfill(3)))
        plt.close()

        plt.scatter(cosmos_rmag-cosmos_imag, 
                    cosmos_gmag-cosmos_imag, 
                    edgecolor='none', c='k', label='all cosmos high z')
        plt.scatter(cosmos_rmag[inear]-cosmos_imag[inear], 
                    cosmos_gmag[inear]-cosmos_rmag[inear], 
                    edgecolor='none', c='b', label='15 nearest cosmos high z')
        plt.scatter(-2.5*np.log10(r[wrong][icheck]/i[wrong][icheck]),
                    -2.5*np.log10(g[wrong][icheck]/r[wrong][icheck]),
                    edgecolor='none', c='m')
        plt.xlabel('r - i')
        plt.ylabel('g - r')
        plt.title('truth z = {}, Phigh = {}'.format(truthz[wrong][icheck], Phigh[wrong][icheck]))
        plt.xlim(-1,4)
        plt.ylim(-1,4)
        plt.legend(loc='best')
        plt.savefig('{}_colors_imposter_{}'.format(base_name, str(icheck).zfill(3)))
        plt.close()

    for ni in range(ncheck):
        plt.plot(nfilt, [g[wrong][ni], r[wrong][ni], i[wrong][ni], z[wrong][ni]], c='r')
        plt.plot(nfilt, [g[right][ni], r[right][ni], i[right][ni], z[right][ni]], c='g')
    plt.xlabel('filter')
    plt.ylabel('flux_detmodel')
    plt.title('red=wrong, green=right, cutoff: P>'+str(cut))
    plt.savefig('{}_wrongright_fluxes'.format(base_name))
    plt.close()
    
#plot fluxes of nearest COSMOS

if __name__=="__main__":
    main()

    
