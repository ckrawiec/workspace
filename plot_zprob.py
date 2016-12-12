import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, join
from astropy.io import fits
import glob
import pandas as pd
from myutils import fitsstack

base_name = 'zprob_SV_z25_3bins_sigma_tree_griz'
#'zprob_balrog_sva1_z25_3bin_sigma_tree_griz'
#'zprob_Y1_z25_3bins_sigma_tree_griz'
#'zprob_balrog_y1a1_z25_3bins_sigma_tree_griz'


zsrc = [2.5, 9.9]
zlens = [0.001, 0.8]
zother = [0.8,2.5]

d_file = '/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_good_regions_no_cosmos.fits'
#'/home/ckrawiec/DES/data/y1a1_gold_mag_detmodel_MC1.fits'
#'/home/ckrawiec/DES/data/balrog_y1a1_truth_sim_flux_detmodel.fits'

red_file = '/home/ckrawiec/DES/data/y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10.fit'

z_file = '/home/ckrawiec/DES/magnification/lbgselect/{}.fits'.format(base_name)
z_files = glob.glob('/home/ckrawiec/DES/magnification/lbgselect/{}_*.fits'.format(base_name))

def makeplots():
#    hexbinPsrcmean()
#    hist2dPall()
#    NvsPcut()
#    radecscatter(0.6)
#    histmagi(0.6)
    histmagislope(0.6)
#    redmagiccount(0.6)
#    histhexbin()
#    hexbinPsrcmeanNorm()
#    balrogzhist(0.8)
        
def writefitsstack(infiles, outfile):
    hdu = fitsstack(infiles)
    pri_hdu = fits.PrimaryHDU()
    hdu_list = fits.HDUList([pri_hdu, hdu])
    hdu_list.writeto(outfile)

#writefitsstack(z_files, z_file)

ztab = Table.read(z_file)
dtab = Table.read(d_file)
#zbalrogtab = Table.read('/home/ckrawiec/DES/data/balrog_y1a1_truth_index_z.fits')

#newtab = join(dtab, zbalrogtab)
tab = join(dtab, ztab)

del ztab, dtab#, zbalrogtab

if 'balrog' not in base_name:
    g = tab['MAG_DETMODEL_G']
    r = tab['MAG_DETMODEL_R']
    i = tab['MAG_DETMODEL_I']
    z = tab['MAG_DETMODEL_Z']
    Y = tab['MAG_DETMODEL_Y']
    
    grbox = (((g-r) < 4.) & ((g-r) > -1.5))
    ribox = (((r-i) < 4.) & ((r-i) > -1.5))
    izbox = (((i-z) < 4.) & ((i-z) > -1.5))
    zYbox = (((z-Y) < 4.) & ((z-Y) > -1.5))
    
    gr = g-r
    ri = r-i
    iz = i-z
    zY = z-Y

    del g,r,z,Y
    
ra, dec = tab['RA'], tab['DEC']

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

def histmagislope(cut):
    h = plt.hist(i[spte & (Psrc>cut) & (i>17) & (i<28)], histtype='step', bins=20, label='SPT-E')
    plt.yscale('log')

    dx = np.gradient(h[1][:-1])
    logs = np.log10(h[0])
    alphas = 2.5 * np.gradient(logs, dx)
    
    plt.plot(h[1][:-1], alphas, label='alpha')
    plt.xlim(17,28)
    plt.legend(loc='upper left')
    plt.xlabel('mag_detmodel_i')
    plt.title('P(z > {}) > {}'.format(zsrc[0], cut))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_hist_mag_i_diffslope_Pcut_src'.format(base_name))
    plt.close()


def histmagi(cut):
    plt.hist(i[Psrc>cut], histtype='step', bins=1000, label='all')
    plt.hist(i[spte & (Psrc>cut)], histtype='step', bins=1000, label='SPT-E')
#    plt.hist(i[cosmos & (Psrc>cut)], histtype='step', bins=1000, label='COSMOS')
#    plt.hist(i[y1main & (Psrc>cut)], histtype='step', bins=1000, label='y1a1 main')
    plt.yscale('log')
    plt.xlim(17,26)
    plt.legend(loc='upper left')
    plt.xlabel('mag_detmodel_i')
    plt.title('P(z > {}) > {}'.format(zsrc[0], cut))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_hist_mag_i_Pcut_src'.format(base_name))
    plt.close()

    plt.hist(i[Plens>cut], histtype='step', bins=1000, label='all')
    plt.hist(i[spte & (Plens>cut)], histtype='step', bins=1000, label='SPT-E')
#    plt.hist(i[cosmos & (Psrc>cut)], histtype='step', bins=1000, label='COSMOS')
    plt.hist(i[y1main & (Plens>cut)], histtype='step', bins=1000, label='y1a1 main')
    plt.yscale('log')
    plt.xlim(17,26)
    plt.legend(loc='upper left')
    plt.xlabel('mag_detmodel_i')
    plt.title('P({} < z < {}) > {}'.format(zlens[0], zlens[1], cut))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_hist_mag_i_Pcut_lens'.format(base_name))
    plt.close()


def histhexbin():
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df.plot.hexbin(x='r-i', y='g-r', gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_ri_gr_hist_hexbin'.format(base_name))
    plt.close()

    df = pd.DataFrame(zip(iz[izbox & ribox], ri[izbox & ribox]), columns=['i-z', 'r-i'])
    df.plot.hexbin(x='i-z', y='r-i', gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_iz_ri_hist_hexbin'.format(base_name))
    plt.close()
    
    df = pd.DataFrame(zip(zY[zYbox & izbox], iz[zYbox & izbox]), columns=['z-Y', 'i-z'])
    df.plot.hexbin(x='z-Y', y='i-z', gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_zY_iz_hist_hexbin'.format(base_name))
    plt.close()

def colorfunc(color):
    return np.log10(np.mean(color)*len(color))
    
def hexbinPsrcmeanNorm():
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df['Psrc']=Psrc[ribox & grbox]
    df.plot.hexbin(x='r-i', y='g-r', C='Psrc', reduce_C_function=colorfunc, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_ri_gr_P_mean_norm_hexbin'.format(base_name))
#    plt.colorbar(label='P{} * N'.format(zsrc))
    plt.close()

def hexbinPsrcmean(): 
    
    df = pd.DataFrame(zip(ri[ribox & grbox], gr[ribox & grbox]), columns=['r-i', 'g-r'])
    df['Psrc']=Psrc[ribox & grbox]
    df.plot.hexbin(x='r-i', y='g-r', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_ri_gr_P_mean_hexbin'.format(base_name))
    plt.close()

    df = pd.DataFrame(zip(iz[izbox & ribox], ri[izbox & ribox]), columns=['i-z', 'r-i'])
    df['Psrc']=Psrc[izbox & ribox]
    df.plot.hexbin(x='i-z', y='r-i', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_iz_ri_P_mean_hexbin'.format(base_name))
    plt.close()
    
    df = pd.DataFrame(zip(zY[zYbox & izbox], iz[zYbox & izbox]), columns=['z-Y', 'i-z'])
    df['Psrc']=Psrc[zYbox & izbox]
    df.plot.hexbin(x='z-Y', y='i-z', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('gnuplot'))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_zY_iz_P_mean_hexbin'.format(base_name))
    plt.close()

    
def hist2dPall():
    plt.hist2d(Plens[notnan], Psrc[notnan], bins=100, norm=mpl.colors.LogNorm())
    plt.colorbar()
    plt.xlabel('P'+str(zlens))
    plt.ylabel('P'+str(zsrc))
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_Pplot'.format(base_name))
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
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_N_Pcut'.format(base_name))
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
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_ra_dec'.format(base_name))

    #Y1A1 main
    plt.xlim(0,360)
    plt.ylim(-70,-35)
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_ra_dec_Y1A1_main'.format(base_name))

    #SPT-E
    plt.xlim(60, 95)
    plt.ylim(-62.3, -42)
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_ra_dec_SPTE'.format(base_name))

    #COSMOS
    plt.xlim(148.5, 151.5)
    plt.ylim(1., 3.5)
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_ra_dec_COSMOS'.format(base_name))
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
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_redmagic_annulus_count'.format(base_name))
    plt.close()

def balrogzhist(cut):
    plt.hist(redshift[Psrc>cut], histtype='step',label='P{} > {}'.format(zsrc, cut), bins=100)
    plt.hist(redshift[Plens>cut], histtype='step',label='P{} > {}'.format(zlens, cut), bins=100)
    plt.hist(redshift[Pother>cut], histtype='step',label='P{} > {}'.format(zother, cut), bins=100)
    plt.yscale('log')
    plt.xlabel('true z')
    plt.legend(loc='best')
    plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_balrog_z_hist'.format(base_name))
    plt.close()


makeplots()
