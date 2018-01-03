import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join, Column
from myutils import writefitsstack, joincorrecttype, printtime

base_name = 'sva1_gold_cosmos_zminchi2_auto_griz_z3_3bins_full_gauss'
#'zprob_SV_z25_3bins_sigma_tree_griz'
#'zprob_balrog_sva1_z25_3bin_sigma_tree_griz'
#'zprob_Y1_z25_3bins_sigma_tree_griz'
#'zprob_balrog_y1a1_z25_3bins_sigma_tree_griz'

home_dir = '/home/ckrawiec/'
SE_col = 'AUTO'

z_groups = [[0.001, 1.0],
            [1.0, 3.0],
            [3.0, 9.9]]
zsrc = [3.0, 9.9]
zlens = [0.001, 1.0]
zother = [1.0, 3.0]

d_id_col = 'COADD_OBJECTS_ID'
z_id_col = 'COADD_OBJECTS_ID'
id_col_type = float

d_file = home_dir+'DES/data/sva1_gold_auto_good_regions_no_cosmos.fits'
#'/home/ckrawiec/DES/data/y1a1_gold_mag_detmodel_MC1.fits'
#'/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel.fits'

red_file = home_dir+'DES/data/y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10.fit'

z_file = home_dir+'DES/magnification/lbgselect/zproboutput/{}_combined.fits'.format(base_name)
z_files = glob.glob(home_dir+'DES/magnification/lbgselect/zproboutput/{}_*.fits'.format(base_name))

def makeplots():
#    hexbinPsrcmean()
#    hist2dPall()
#    NvsPcut()
#    radecscatter(0.6)
#    histmagi(0.9)
    #histmagislope(0.9)
#    redmagiccount(0.6)
#    histhexbin()
    colorscatter()
#    hexbinPsrcmeanNorm()
#    balrogzhist(0.8)

printtime(stderr=True)

if len(z_files) == 0:
    sys.stderr.write("chosen base name resulted in no results. exiting.\n")
    exit()

if os.path.exists(z_file):
    print "combined zprob file already exists! ({})".format(z_file)
    print "continuing..."
else:
    writefitsstack(z_files, z_file)
    print "successfully combined zprob files"

#ztab = Table.read(z_file)
#dtab = Table.read(d_file)
#zbalrogtab = Table.read('/home/ckrawiec/DES/data/balrog_y1a1_truth_index_z.fits')
#newtab = join(dtab, zbalrogtab)
#del ztab, dtab#, zbalrogtab
printtime(stderr=True)
sys.stderr.write("before join\n")
tab = joincorrecttype(d_file, z_file, d_id_col, z_id_col, id_col_type)
printtime(stderr=True)
sys.stderr.write("after join\n")

if 'MAG_'+SE_col+'_G' in tab.colnames:
    g = tab['MAG_'+SE_col+'_G']
    r = tab['MAG_'+SE_col+'_R']
    i = tab['MAG_'+SE_col+'_I']
    z = tab['MAG_'+SE_col+'_Z']
    Y = tab['MAG_'+SE_col+'_Y']

    gr = g-r
    ri = r-i
    iz = i-z
    zY = z-Y
    
    del g,r,z,Y
    
    grbox = (((gr) < 4.) & ((gr) > -1.5))
    ribox = (((ri) < 4.) & ((ri) > -1.5))
    izbox = (((iz) < 4.) & ((iz) > -1.5))
    zYbox = (((zY) < 4.) & ((zY) > -1.5))

    
ra, dec = tab['RA'], tab['DEC']

P = {}
for z_group in z_groups:
    P[str(z_group)] = tab['P'+str(z_group)]
Psrc = tab['P'+str(zsrc)]
Plens = tab['P'+str(zlens)]
#Pother = tab['P'+str(zother)]

if 'balrog' in base_name:
    redshift = tab['Z']
    
del tab

cosmos = (ra>148.5) & (ra>151.5) & (dec>1.) & (dec<3.5)
spte = (ra>50) & (ra<99) & (dec>-65) & (dec<-40)
y1main = (dec < -35)

notnan = ~np.isnan(Plens) & ~np.isnan(Psrc)

def colorscatter():
    printtime(stderr=True)
    plt.scatter(ri[ribox & grbox], gr[ribox & grbox], c=Psrc[ribox & grbox], edgecolor='none', s=2.)
    plt.xlabel('r - i')
    plt.ylabel('g - r')
    plt.colorbar(label='P(z in {})'.format(zsrc))
    plt.xlim(-1.5, 4.)
    plt.ylim(-1.5, 4.)
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_Psrc_colorscatter.png'.format(base_name))
    plt.close()

    plt.scatter(ri[ribox & grbox], gr[ribox & grbox], c=Plens[ribox & grbox], edgecolor='none', s=2.)
    plt.xlabel('r - i')
    plt.ylabel('g - r')
    plt.colorbar(label='P(z in {})'.format(zlens))
    plt.xlim(-1.5, 4.)
    plt.ylim(-1.5, 4.)
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_Plens_colorscatter.png'.format(base_name))
    plt.close()

def histmagislope(cut):
    printtime(stderr=True)
    sys.stderr.write("working on i-magnitude histogram slopes...\n")
    for z_group in z_groups:
        h = plt.hist(i[spte & (P[str(z_group)]>cut) & (i>17) & (i<30)], 
                     histtype='step', bins=26, label='SPT-E')
        if len(h[1][:-1])>0:
            plt.yscale('log')
        
            dx = np.gradient(h[1][:-1])
            logs = np.log10(h[0])
            alphas = 2.5 * np.gradient(logs, dx)

            plt.plot(h[1][:-1], alphas, label='alpha')
            plt.xlim(17,28)
            plt.legend(loc='upper left')
            plt.xlabel('mag_detmodel_i')
            plt.title('P(z in {}) > {}'.format(z_group, cut))
            plt.savefig(home_dir+'DES/magnification/lbgselect/{}_hist_mag_i_diffslope_Pcut{}_{}.png'.format(base_name, 
                                                                                                            cut, 
                                                                                                            z_group))
        plt.close()

    sys.stderr.write("done\n")


def histmagi(cut):
    printtime(stderr=True)
    sys.stderr.write("working on i-magnitude histograms...\n")
    nbins = 26
    imask = ((i>17) & (i<30))

    for z_group in z_groups:
        plt.hist(i[(P[str(z_group)]>cut) & imask & spte], histtype='step', bins=nbins, normed=True,
                 label='P(z in {}) > {}'.format(z_group, cut))
#        plt.hist(i[spte & (P[str(z_group)]>cut) & imask], histtype='step', bins=nbins, label='SPT-E')
#        plt.hist(i[cosmos & (P[str(z_group)]>cut) & imask], histtype='step', bins=nbins, label='COSMOS')
#        plt.hist(i[y1main & (P[str(z_group)]>cut)], histtype='step', bins=nbins, label='y1a1 main')
    plt.yscale('log')
    plt.xlim(17,30)
    plt.legend(loc='upper left')
    plt.xlabel('MAG_'+SE_col+'_I')
    plt.ylabel('normalized counts')
    plt.title('SPT-E')
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_hist_mag_i_Pcut{}.png'.format(base_name, 
                                                                                           cut))
    plt.close()

    sys.stderr.write("done\n")


def histhexbin():
    sys.stderr.write("working on hexbin histograms...\n")
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df.plot.hexbin(x='r-i', y='g-r', gridsize=20, bins='log', cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_ri_gr_hist_hexbin'.format(base_name))
    plt.close()

    df = pd.DataFrame(zip(iz[izbox & ribox], ri[izbox & ribox]), columns=['i-z', 'r-i'])
    df.plot.hexbin(x='i-z', y='r-i', gridsize=20, bins='log', cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_iz_ri_hist_hexbin'.format(base_name))
    plt.close()
    
    df = pd.DataFrame(zip(zY[zYbox & izbox], iz[zYbox & izbox]), columns=['z-Y', 'i-z'])
    df.plot.hexbin(x='z-Y', y='i-z', gridsize=20, bins='log', cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_zY_iz_hist_hexbin'.format(base_name))
    plt.close()
    sys.stderr.write("done\n")

def colorfunc(color):
    return np.log10(np.mean(color)*len(color))
    
def hexbinPsrcmeanNorm():
    sys.stderr.write("working on hexbin means...\n")
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df['Psrc']=Psrc[ribox & grbox]
    df.plot.hexbin(x='r-i', y='g-r', C='Psrc', reduce_C_function=colorfunc, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_ri_gr_Psrc_mean_norm_hexbin'.format(base_name))
#    plt.colorbar(label='P{} * N'.format(zsrc))
    plt.close()

def hexbinPsrcmean(): 
    
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df['Psrc']=Psrc[ribox & grbox]
    df.plot.hexbin(x='r-i', y='g-r', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_ri_gr_P_mean_hexbin'.format(base_name))
    plt.close()

    df = pd.DataFrame(zip(iz[izbox & ribox], ri[izbox & ribox]), columns=['i-z', 'r-i'])
    df['Psrc']=Psrc[izbox & ribox]
    df.plot.hexbin(x='i-z', y='r-i', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_iz_ri_P_mean_hexbin'.format(base_name))
    plt.close()
    
    df = pd.DataFrame(zip(zY[zYbox & izbox], iz[zYbox & izbox]), columns=['z-Y', 'i-z'])
    df['Psrc']=Psrc[zYbox & izbox]
    df.plot.hexbin(x='z-Y', y='i-z', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_mag_zY_iz_P_mean_hexbin'.format(base_name))
    plt.close()
    sys.stderr.write("done\n")
    
def hist2dPall():
    plt.hist2d(Plens[notnan], Psrc[notnan], bins=100, norm=mpl.colors.LogNorm())
    plt.colorbar()
    plt.xlabel('P'+str(zlens))
    plt.ylabel('P'+str(zsrc))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_Pplot'.format(base_name))
    plt.close()

def NvsPcut():
    Pcut = np.arange(0.1, 1.0, 0.01)
    nsrc = np.array([len(Psrc[Psrc>icut]) for icut in Pcut])
    nlens = np.array([len(Plens[Plens<icut]) for icut in Pcut])
    ncosmos = np.array([len(Psrc[cosmos & (Psrc>icut)]) for icut in Pcut])
    nspte = np.array([len(Psrc[spte & (Psrc>icut)]) for icut in Pcut])
    ny1main = np.array([len(Psrc[y1main & (Psrc>icut)]) for icut in Pcut])
    
#    plt.plot(Pcut, ncosmos, c='g', label='$P(z>'+str(zsrc[0])+') > P_{cut}$ in COSMOS area')
    plt.plot(Pcut, nspte, c='r', label='$P(z>'+str(zsrc[0])+') > P_{cut}$ in SPT-E area')
#    plt.plot(Pcut, ny1main, c='b', label='$P(z>'+str(zsrc[0])+') > P_{cut}$ in Y1A1 main area')
#    plt.plot(Pcut, nsrc, c='k', label='$P(z>'+str(zsrc[0])+') > P_{cut}$')
    plt.yscale('log')
    plt.xlabel('$P_{cut}$')
    plt.ylabel('N')
    plt.legend(loc='best')
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_N_Pcut'.format(base_name))
    plt.close()
    
def radecscatter(cut):
    hicut = (Psrc > cut)
    locut = (Plens > cut)

    plt.scatter(ra[~hicut], dec[~hicut], edgecolor='none', s=4., label='all')
    plt.scatter(ra[hicut], dec[hicut],c='r', edgecolor='none', s=4., label='P(z>{})>'.format(zsrc[0])+str(cut))
    plt.legend(loc='best')
    plt.xlabel('ra')
    plt.ylabel('dec')

    #all
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_ra_dec'.format(base_name))

    #Y1A1 main
    plt.xlim(0,360)
    plt.ylim(-70,-35)
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_ra_dec_Y1A1_main'.format(base_name))

    #SPT-E
    plt.xlim(60, 95)
    plt.ylim(-62.3, -42)
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_ra_dec_SPTE'.format(base_name))

    #COSMOS
    plt.xlim(148.5, 151.5)
    plt.ylim(1., 3.5)
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_ra_dec_COSMOS'.format(base_name))
    plt.close()

def redmagiccount(cut):
    red_tab = Table.read(red_file)
    redra, reddec = red_tab['RA'], red_tab['DEC']

    r_min = 0.1
    r_max = 1.
    n_src = []
    for this_ra, this_dec in zip(redra,reddec):
        r = np.sqrt((ra-this_ra)**2. + (dec-this_dec)**2.)
        annulus = (r <  r_max) & (r > r_min) & (Psrc > cut)
        n_src.append(len(ra[annulus]))
    plt.hist(n_src, bins=100, histtype='step')
    plt.xlabel('# P(src) > {} around redmagic galaxies ({}-{}$\deg$)'.format(cut,r_min,r_max))
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_redmagic_annulus_count'.format(base_name))
    plt.close()

def balrogzhist(cut):
    plt.hist(redshift[Psrc>cut], histtype='step',label='P{} > {}'.format(zsrc, cut), bins=100)
    plt.hist(redshift[Plens>cut], histtype='step',label='P{} > {}'.format(zlens, cut), bins=100)
    plt.hist(redshift[Pother>cut], histtype='step',label='P{} > {}'.format(zother, cut), bins=100)
    plt.yscale('log')
    plt.xlabel('true z')
    plt.legend(loc='best')
    plt.savefig(home_dir+'DES/magnification/lbgselect/{}_balrog_z_hist'.format(base_name))
    plt.close()


makeplots()
