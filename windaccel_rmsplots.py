#!usr/bin/python
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

accel = ascii.read('reaccelerometeranalysis/accel_rms.txt', guess=False, delimiter=' ')
tab = ascii.read('anemom042114_051514.csv')

accel['speed_rms'] = [0.]*len(accel)
accel['mag_rms'] = [0.]*len(accel)
accel['exptime'] = [0.]*len(accel)

i = 0
j = 0
while j < len(accel):
    if accel['EXPID'][j] == tab['expid'][i]:
        k = i
        while accel['EXPID'][j] == tab['expid'][k]:
            accel['speed_rms'][j] = tab['rms_speed_exp'][k]
            accel['mag_rms'][j] = tab['rms_mag_exp'][k]
            accel['exptime'][j] = tab['exptime'][k]
            k += 1
            if k == len(tab)-1:
                break
            i = k
    else:
        i+=1
        continue
#    print accel['EXPID'][j], ': done!'
    j+=1


jitter = np.sqrt(accel['X-accel(G-rms)']**2 + accel['Y-accel(G-rms)']**2 + accel['Z-accel(G-rms)']**2)
xacc = accel['X-accel(G-rms)']
yacc = accel['Y-accel(G-rms)']
zacc = accel['Z-accel(G-rms)']

"""
plt.scatter(accel['speed_rms'], accel['mag_rms'])
plt.xlabel('rms speed')
plt.ylabel('rms magnitude')
plt.savefig('magrms_vs_speedrms')
plt.close()

plt.scatter(accel['speed_rms'], jitter, c='r')
plt.ylabel('rms total acceleration (G)')
plt.xlabel('rms wind speed')
plt.ylim(2.8e-4,4e-4)
plt.savefig('totaccel_vs_speedrms')
plt.close()
"""

varsub = jitter**2 - np.median(jitter**2)


"""
for n in range(len(accel)):
   if accel['exptime'][n] == 15.:
        plt.scatter(accel['mag_rms'][n], jitter[n]**2 - np.median(jitter**2), c='orange', edgecolor='orange')
   else:
        plt.scatter(accel['mag_rms'][n], jitter[n]**2 - np.median(jitter**2), c='g', edgecolor='g')
plt.ylim(-.2e-7,.6e-7)
plt.xlabel('rms magnitude')
plt.ylabel('median-subtracted mean square tot. acceleration by exposure time')
plt.savefig('excessaccel_vs_magrms_byexp')
plt.close()

ptot = plt.scatter(accel['speed_rms'], jitter**2 - np.median(jitter**2))
px = plt.scatter(accel['speed_rms'], xacc**2 - np.median(xacc**2), c='r')
py = plt.scatter(accel['speed_rms'], yacc**2 - np.median(yacc**2), c='y')
pz = plt.scatter(accel['speed_rms'], zacc**2 - np.median(zacc**2), c = 'g')
x = np.linspace(0, 9)
y = [1.e-8]*len(x)
plt.plot(x,y, c='black')
plt.legend([ptot, px, py, pz], ['Total-accel','X-accel', 'Y-accel', 'Z-accel'], loc='upper left')
plt.xlabel('rms w_speed')
plt.ylabel('median-subtracted acceleration variance (G^2)')
plt.ylim(-.2e-7,.6e-7)
plt.xlim(0,4)
plt.savefig('excessaccel_vs_speedrms')
plt.close()
"""

ptot = plt.scatter(accel['mag_rms'], jitter**2 - np.median(jitter**2))
px = plt.scatter(accel['mag_rms'], xacc**2 - np.median(xacc**2), c='r', edgecolor='none')
py = plt.scatter(accel['mag_rms'], yacc**2 - np.median(yacc**2), c='y', edgecolor='none')
pz = plt.scatter(accel['mag_rms'], zacc**2 - np.median(zacc**2), c = 'g', edgecolor='none')
plt.legend([ptot, px, py, pz], ['Total-accel','X-accel', 'Y-accel', 'Z-accel'], loc='upper left')
x = np.linspace(0, 9)
y = [1.e-8]*len(x)
plt.plot(x,y,c='black')
plt.text(2, y[0], '$10^{-8} G^{2}$')
plt.xlabel('rms wind speed (m/s)') #rms magnitude
plt.ylabel('median-subtracted acceleration variance ($G^{2}$)')
plt.ylim(-.2e-7,.6e-7)
plt.xlim(0,5)
plt.savefig('excessaccel_vs_magrms2')
plt.close()

"""
plt.scatter(accel['mag_rms'], jitter, c='g')
#plt.plot(np.polyfit(accel['mag_rms'], jitter, 1))
plt.ylabel('rms total acceleration (G)')
plt.xlabel('rms magnitude')
plt.ylim(1e-4,4e-4)
plt.savefig('totaccel_vs_magrms')
plt.close()
"""


#p11 = plt.scatter(accel['speed_rms'], jitter, c='r')
#p22 = plt.scatter(accel['mag_rms'], jitter, c='g')
p33 = plt.scatter(np.sqrt(accel['mag_rms']**2+accel['speed_rms']**2), jitter, edgecolor='none')
#plt.legend([p11,p22,p33], ['rms speed', 'rms magnitude', 'rms total'])
plt.ylabel('total rms acceleration (G)')
plt.xlabel('rms wind speed (m/s)')
plt.ylim(2.9e-4,4e-4)
plt.savefig('totaccel_vs_windspeeds2')
plt.close()


plt.scatter(np.sqrt(accel['mag_rms']**2+accel['speed_rms']**2), jitter)
plt.ylabel('rms total acceleration (G)')
plt.xlabel('wind')
plt.ylim(2.8e-4,4e-4)
plt.savefig('totaccel_vs_totwspeed2')
plt.close()

p1 = plt.scatter(accel['mag_rms'], accel['X-accel(G-rms)'])
p2 = plt.scatter(accel['mag_rms'], accel['Y-accel(G-rms)'], c='g')
p3 = plt.scatter(accel['mag_rms'], accel['Z-accel(G-rms)'], c='r')
plt.legend([p1, p2, p3], ['X-accel', 'Y-accel', 'Z-accel'])
plt.ylabel('accel')
plt.xlabel('rms magnitude')
plt.ylim(1e-4,3e-4)
plt.savefig('acceldirs_vs_magrms2')
plt.close()
"""

jwjoin = zip(varsub, accel['mag_rms'], accel['EXPID'])

f = open('accel_outliers','w')

for item in jwjoin:
    if item[0] >= 1.e-8 and item[1] < 2.:
        f.write(str(item[2])+'\n')

f.close()


plt.plot(accel['EXPID'], accel['mag_rms'], 'go-')
plt.xlabel('exposure id')
plt.ylabel('rms magnitude (m/s)')
plt.savefig('magrms_vs_exp')
plt.close()

plt.plot(accel['EXPID'], jitter, 'bo-') 
plt.xlabel('exposure id')
plt.ylabel('total rms acceleration (G)')
plt.ylim(2.9e-4,4e-4)
plt.savefig('rmsaccel_vs_exp')
plt.close()

exps = [item['expid'] for item in tab if item['expid']!=0]
dates = [item['time_recorded'][6:10] for item in tab if item['expid']!=0]
mags = [item['rms_mag_exp'] for item in tab if item['expid']!=0]

x = range(0,len(dates))
tix = zip(x,dates)
timetix = [item[1] for item in tix if item[0] in np.arange(0,len(exps),15000)]
plt.xticks(np.arange(0,len(exps),15000), timetix, rotation=30)
plt.plot(x, mags, 'yo-')
#plt.plot(exps, mags, 'yo-')
plt.xlabel('date')
plt.ylabel('rms magnitude')
plt.savefig('rmsmag_vs_date_421_515')
plt.close()


host = host_subplot(111, axes_class=AA.Axes)
plt.subplots_adjust(right=0.75)
par1 = host.twinx()
host.set_xlabel("Exposure ID")
host.set_ylabel("Total Rms Acceleration (G)")
par1.set_ylabel("Rms Magnitude (m/s)")
p2, = host.plot(accel['EXPID'],jitter, 'bo-')
p1, = par1.plot(accel['EXPID'],accel['mag_rms'], 'ro-')
host.set_ylim(2.99e-4)
par1.set_ylim(0., 5.)
host.axis["left"].label.set_color(p2.get_color())
par1.axis["right"].label.set_color(p1.get_color())
plt.draw()
plt.savefig('rmsaccelmag_vs_exp')
"""
exit()


