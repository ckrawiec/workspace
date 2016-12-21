import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, join

zprob_file = '/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_z25_3bins_sigma_tree_griz_tab{}.fits'
sim_file = '/home/ckrawiec/DES/data/balrog_sva1_tab{}_SIM.fits'
truth_file = '/home/ckrawiec/DES/data/balrog_sva1_tab{}_TRUTH.fits'

num_list = [str(i).zfill(2) for i in range(1,12)]
num_list.remove('09')

zsrc = [2.5, 9.9]
zlens = [0.001, 0.8]
zother = [0.8, 2.5]

Plow, Phigh, Pother, truthz, ra, dec = [],[],[],[],[],[]
g,r,i,z = [],[],[],[]

def main():
    #plotStats()
    printStats()
    #radec()
    #Phightruthz()
    #truthzhist(0.6)
    checkImposters(cut=0.8)

for num in num_list:
    zprob = Table.read(zprob_file.format(num))
    sim = Table.read(sim_file.format(num))
    truth = Table.read(truth_file.format(num))

    tab1 = join(zprob, sim)
    tab = join(tab1, truth)

    Plow.append(tab['P'+str(zlens)])
    Phigh.append(tab['P'+str(zsrc)])
    Pother.append(tab['P'+str(zother)])

    truthz.append(tab['Z'])

    ra.append(tab['RA'])
    dec.append(tab['DEC'])

    g.append(tab['FLUX_DETMODEL_G'])
    r.append(tab['FLUX_DETMODEL_R'])
    i.append(tab['FLUX_DETMODEL_I'])
    z.append(tab['FLUX_DETMODEL_Z'])

g = np.hstack(g)
r = np.hstack(r)
i = np.hstack(i)
z = np.hstack(z)

Plow = np.hstack(Plow)
Phigh = np.hstack(Phigh)
Pother = np.hstack(Pother)
truthz = np.hstack(truthz)

def plotStats():
    truthzCutHist(Plow, zlens)
    truthzCutHist(Phigh, zsrc)
    truthzCutHist(Pother, zother)

def truthzCutHist(Pgrp, zgrp):
    zmask = (truthz>np.min(zgrp)) & (truthz<np.max(zgrp))
    for cut in np.arange(0,1.0,0.2):
        plt.hist(truthz[zmask & (Pgrp>cut)], histtype='step', label='P({})>{}'.format(zgrp, cut))
    plt.xlabel('truth z')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/balrog_sva1_truthz_{}cut_hist.png'.format(zgrp))
    plt.close()

    for cut in np.arange(0,1.0,0.2):
        plt.hist(truthz[Pgrp>cut], histtype='step', label='P({})>{}'.format(zgrp, cut))
    plt.xlabel('truth z')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/balrog_sva1_truthz_{}cut_hist_allz.png'.format(zgrp))
    plt.close()

def printStats():
    print "Number of galaxies with truth z > 2.5: ", len(truthz[truthz>2.5])
    print "Number of galaxies with truth z > 2.5 & Phigh > 0.7: ", len(truthz[(truthz>2.5) & (Phigh>0.7)])
    print "Number of galaxies with Phigh > 0.7: ", len(truthz[Phigh>0.7])

def radec():    
    plt.scatter(np.hstack(ra), np.hstack(dec), edgecolor='none')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.savefig('/home/ckrawiec/DES/data/balrog_sva1_SIM_ra_dec')
    plt.close()
    
def Phightruthz():
    plt.scatter(truthz, Phigh,
                s=4., edgecolor='none', c=Plow)
    plt.colorbar(label='P[0.001, 0.8]')
    plt.xlabel('truth z')
    plt.ylabel('P[2.5, 9.9]')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_z25_3bins_sigma_tree_griz_truthz_Phigh_Plow_scatter')
    plt.close()

    notnan = (~np.isnan(Plow)) & (~np.isnan(Phigh)) & (~np.isnan(Pother))

    plt.hist2d(truthz[notnan], Phigh[notnan], bins=100.)
    plt.xlabel('truth z')
    plt.ylabel('P[2.5, 9.9]')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_z25_3bins_sigma_tree_griz_truthz_Phigh_hist2d')
    plt.close()

def truthzhist(cut):
    plt.hist(truthz[Phigh>cut], histtype='step',label='P{} > {}'.format(zsrc, cut), bins=100)
    plt.hist(truthz[Plow>cut], histtype='step',label='P{} > {}'.format(zlens, cut), bins=100)
    plt.hist(truthz[Pother>cut], histtype='step',label='P{} > {}'.format(zother, cut), bins=100)
    plt.yscale('log')
    plt.xlabel('true z')
    plt.legend(loc='best')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_z25_3bins_sigma_tree_griz_z_hist')
    plt.close()

def checkImposters(cut):
    wrong = (Phigh>cut) & (truthz<zsrc.min()) & (truthz>zsrc.max())
    right = (Phigh>cut) & (truthz>zsrc.min()) & (truthz>zsrc.min())
    nfilt = range(4)
    ncheck = 10
    for ni in range(ncheck):
        plt.plot(nfilt, [g[wrong][ni], r[wrong][ni], i[wrong][ni], z[wrong][ni]], c='r')
        plt.plot(nfilt, [g[right][ni], r[right][ni], i[right][ni], z[right][ni]], c='g')
    plt.xlabel('filter')
    plt.ylabel('flux_detmodel')
    plt.title('red=wrong, green=right, cutoff: P>'+str(cut))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/zprob_balrog_sva1_z25_3bins_sigma_tree_griz_wrongright_fluxes')
    plt.close()
    
#plot fluxes of nearest COSMOS

if __name__=="__main__":
    main()

    
