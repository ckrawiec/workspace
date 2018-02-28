import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import glob

tab = ascii.read('complete_trunc_rms_stats.tab')
tel9 = ascii.read('tel_091714_101714_updated.csv')

total_bins = 10
mag = tab['rms_mag_exp']
bins = np.linspace(mag.min(),mag.max(), total_bins)
df = pd.DataFrame({'X': mag, 'Y': tab['atotrms'], 'Z': tab['dtotrms'], 'E':tab['exptime']})
data_cut = pd.cut(df.X,bins)
grp = df.groupby(by = data_cut)
ret = grp.aggregate(np.median)
exptimes = tab['exptime']
rain = plt.get_cmap("gist_rainbow")
jet = plt.get_cmap("jet")


plt.scatter(df.X,df.Y,c=exptimes,cmap=rain,edgecolor='none')
plt.plot(ret.X,ret.Y,'ko-',ms=5.)
plt.xlabel('rms wind speed (m/s)')
plt.ylabel('total rms acceleration (m/s^2)')
cbar = plt.colorbar()
cbar.set_clim(0,200)
cbar.set_label('exptime',rotation=270,labelpad=20)
plt.ylim(0.0001,0.0006)
plt.xlim(0,8.1)
fig = plt.gcf()
fig.set_size_inches(9,5)
plt.savefig('rms_totaccel_wind_med2')
plt.close()

plt.scatter(df.X,df.Z*1e6,c=np.log10(df.E),cmap=jet,edgecolor='none')
plt.plot(ret.X,ret.Z*1e6,'ko-',ms=5.)
plt.xlabel('rms wind speed ($m/s$)')
plt.ylabel('total rms displacement ($\mu m$)')
cbar = plt.colorbar()
cbar.set_clim(np.log10(np.min(df.E)),np.log10(np.max(df.E)))
cbar.set_label('log(exposure time ($s$))',rotation=270,labelpad=20)
plt.ylim(0.5,2)
plt.xlim(0,8.1)
fig = plt.gcf()
fig.set_size_inches(10,5)
fig.tight_layout()
plt.savefig('rms_totdisp_wind_med2')
plt.close()

non0 = list(set(np.hstack([np.where(tel9['gyy']!=0.)[0],np.where(tel9['rms_mag_exp']!=0.)[0]])))

magt = tel9['rms_mag_exp'][non0]
dft = pd.DataFrame({'X': tel9['rms_mag_exp'][non0],'Y':tel9['gyy'][non0]})
binst = np.linspace(tel9['rms_mag_exp'][non0].min(),tel9['rms_mag_exp'][non0].max(), total_bins)
cutt = pd.cut(dft.X,binst)
grpt = dft.groupby(by = cutt)
rett = grpt.aggregate(np.median)

gyy = tel9['gyy'][non0]
gxx = tel9['gxx'][non0]
gxy = tel9['gxy'][non0]

pgyy = plt.scatter(magt,gyy,c='r',edgecolor='none')
pyy, = plt.plot(rett.X,rett.Y,'mo-',ms=8)
dft.Y = gxx
grpt = dft.groupby(by = cutt)
rett = grpt.aggregate(np.median)
pxx, = plt.plot(rett.X,rett.Y,'co-',ms=8)
dft.Y = gxy
grpt = dft.groupby(by =cutt)
rett = grpt.aggregate(np.median)
pxy, =plt.plot(rett.X,rett.Y,'yo-',ms=8)
plt.xlabel('rms wind speed (m/s)')
plt.ylabel('guider variance (pixel^2)')
plt.xlim(0,8.)
plt.ylim(-0.005,0.065)
plt.legend([pgyy,pyy,pxx,pxy],['gyy','gyy median','gxx median','gxy median'],loc='upper left')
plt.savefig('guider_med2')
plt.close()

non0 = list(set(np.hstack([np.where(tel9['yy']!=0.)[0],np.where(tel9['rms_mag_exp']!=0.)[0]])))

magt = tel9['rms_mag_exp'][non0]
dft = pd.DataFrame({'X': tel9['rms_mag_exp'][non0],'Y':tel9['yy'][non0]})
binst = np.linspace(tel9['rms_mag_exp'][non0].min(),tel9['rms_mag_exp'][non0].max(), total_bins)
cutt = pd.cut(dft.X,binst)
grpt = dft.groupby(by = cutt)
rett = grpt.aggregate(np.median)

yy = tel9['yy'][non0]
xx = tel9['xx'][non0]
xy = tel9['xy'][non0]

pgxx = plt.scatter(magt,xx,c='b',edgecolor='none')
pyy, =plt.plot(rett.X,rett.Y,'mo-',ms=8)
dft.Y = xx
grpt = dft.groupby(by = cutt)
rett = grpt.aggregate(np.median)
pxx, =plt.plot(rett.X,rett.Y,'co-',ms=8)
dft.Y =xy
grpt = dft.groupby(by =cutt)
rett = grpt.aggregate(np.median)
pxy, =plt.plot(rett.X,rett.Y,'yo-',ms=8)
plt.xlabel('rms wind speed (m/s)')
plt.ylabel('tcs variance (arcsec^2)')
plt.xlim(0,8.)
plt.ylim(-0.025,0.4)
plt.legend([pgxx,pyy,pxx,pxy],['xx','yy median','xx median','xy median'],loc='upper left')
plt.savefig('tcs_med2')
plt.close()

"""
for filename in glob.glob('09172014_062547_PM_subset/unfiltered_accel*fits'):
    dat1 = pf.open(filename)[1]
    dat = dat1.data
    x,y,z = dat['x'],dat['y'],dat['z']
    time = np.arange(0.,len(x)/50.,0.02)
    plt.plot(time,x)
#    plt.plot(time,y)
#    plt.plot(time,z)
    plt.xlabel('Time (s)')
    plt.ylabel('unfiltered voltage')
    plt.title(filename[-11:-5])
    plt.savefig('09172014_062547_PM_subset/unfiltered_'+str(filename[-11:-5]))
    plt.close()
"""
