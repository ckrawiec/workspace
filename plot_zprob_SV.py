import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table, join

base_name = 'zprob_SV_z25_3bins_griz'

zsv = Table.read('/home/ckrawiec/DES/magnification/lbgselect/zprob_SV_new_data_2.fits')
sv = Table.read('/home/ckrawiec/DES/data/sva1_gold_detmodel_MC1_good_regions_no_cosmos.fits')

#zsv.rename_column('coadd_objects_id', 'COADD_OBJECTS_ID')
tab = join(sv, zsv)

del zsv
del sv

zsrc = [2.5, 9.9]
zlens = [0.001, 0.8]
zother = [0.8, 2.5]

g = tab['MAG_DETMODEL_G']
r = tab['MAG_DETMODEL_R']
i = tab['MAG_DETMODEL_I']

Psrc = tab['P'+str(zsrc)]
Plens = tab['P'+str(zlens)]
Pother = tab['P'+str(zother)]

cosmos = (tab['RA']>148.5) & (tab['RA']>151.5) & (tab['DEC']>1.) & (tab['DEC']<3.5)
spte = (tab['RA']>50) & (tab['RA']<99) & (tab['DEC']>-65) & (tab['DEC']<-40)
notnan = ~np.isnan(Plens) & ~np.isnan(Psrc)
box = ((r-i) < 4.) & ((r-i) > -1.) & ((g-r) < 4.) & ((g-r) > -1.)

ri = r[box]-i[box]
gr = g[box]-r[box]

df = pd.DataFrame(zip(ri, gr), columns=['r-i', 'g-r'])
df['Psrc']=Psrc[box]
df.plot.hexbin(x='r-i', y='g-r', C='Psrc', reduce_C_function=np.mean, gridsize=20, cmap=plt.get_cmap('cool'))
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/{}_mag_ri_gr_P_mean_hexbin'.format(base_name))
plt.close()


"""
plt.hist2d(Plens[notnan], Psrc[notnan], bins=100, norm=mpl.colors.LogNorm()) #edgecolor='none', s=4., c=Pother)
plt.colorbar()
plt.xlabel('P'+str(zlens))
plt.ylabel('P'+str(zsrc))
#plt.colorbar(label='P'+str(zother))
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/zprob_SV_z4_2bins_Pplot')
plt.close()
"""

"""                         
Pcut = np.arange(0.1, 1.0, 0.01)
nsrc = np.array([len(Psrc[Psrc>icut]) for icut in Pcut])
nlens = np.array([len(Plens[Plens<icut]) for icut in Pcut])
ncosmos = np.array([len(Psrc[cosmos & (Psrc>icut)]) for icut in Pcut])
nspte = np.array([len(Psrc[spte & (Psrc>icut)]) for icut in Pcut])

plt.plot(Pcut, ncosmos, 'o-', c='g', label='$P(z>'+str(zsrc[0])+') > P_{cut}$ in COSMOS area')
plt.plot(Pcut, nspte, 'o-', c='r', label='$P(z>'+str(zsrc[0])+') > P_{cut}$ in SPT-E area')
plt.plot(Pcut, nsrc, 'o-', c='b', label='$P(z>'+str(zsrc[0])+') > P_{cut}$')
plt.yscale('log')
plt.xlabel('$P_{cut}$')
plt.ylabel('N')
plt.legend(loc='best')
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/N_Pcut_zprob_SV_z4_2bins')
plt.close()
"""

"""
hicut = (Psrc > 0.6)
locut = (Plens > 0.6)

plt.scatter(tab['RA'][~hicut], tab['DEC'][~hicut], edgecolor='none', s=4., label='all')
plt.scatter(tab['RA'][hicut], tab['DEC'][hicut],c='r', edgecolor='none', s=4., label='P(z>{})>0.6'.format(zsrc[0]))
plt.xlim(50, 99)
plt.ylim(-65, -40)
plt.legend(loc='best')
plt.xlabel('ra')
plt.ylabel('dec')
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/ra_dec_zprob_SV_z4_2bins_SPTE')
plt.close()
"""

"""
hicut = (Psrc > 0.8)
locut = (Plens > 0.6)

h0, x0, y0 = np.histogram2d(r[box]-i[box], g[box]-r[box], bins=100)
hhi, xhi, yhi = np.histogram2d(r[box & hicut]-i[box & hicut], g[box & hicut]-r[box & hicut], bins=100)
print 'bins match?: ', np.all(x0==xhi), np.all(y0==yhi)
#hlo = plt.hist2d(r[box & locut]-i[box & locut], g[box & locut]-r[box & locut], bins=h0[1:3])
h = hhi/h0
plt.imshow(h.T, origin='lower', aspect='auto', extent=[-1,4,-1,4], interpolation='nearest', norm=mpl.colors.LogNorm())
plt.xlabel('mag_detmodel_r - mag_detmodel_i')
plt.ylabel('mag_detmodel_g - mag_detmodel_r')
plt.title('bins=100x100')
plt.colorbar(label='# P(z>{})>0.8 / # total'.format(zsrc[0]))
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/mag_ri_gr_hist2d_zprob_SV_z4_2bins')
plt.close()
"""

"""
#plt.scatter(r[box]-i[box], g[box]-r[box], c='k', label='all', edgecolor='none', s=2.)
plt.scatter(r[box & hicut]-i[box & hicut], g[box & hicut]-r[box & hicut],
            c=tab['FLUXERR_DETMODEL_G'][box & hicut], label='P'+str(zsrc)+' > 0.8', edgecolor='none',
            s=4.)
plt.colorbar(label='P'+str(zsrc))
plt.xlim(-1,4)
plt.ylim(-1,4)
plt.xlabel('mag_detmodel_r - mag_detmodel_i')
plt.ylabel('mag_detmodel_g - mag_detmodel_r')
plt.savefig('/home/ckrawiec/DES/magnification/lbgselect/mag_ri_gr_scatter_error_zprob_SV')
plt.close()
"""

#print "# P(source)>0.8: {}, # P(source)>0.8 outside COSMOS field: {}".format(len(Psrc[hicut]), len(Psrc[~cosmos][hicut]))


