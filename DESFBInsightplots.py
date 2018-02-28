import pandas as pd
import pylab as plt
import numpy as np
from astropy.table import Table

f = '/Users/Christina/Documents/DESFBInsightsPosts.csv'
t = '/Users/Christina/Documents/tweet_activity_metrics_theDESurvey_20170726_20170823_en.csv'

#Facebook
reach = 'Lifetime Post Total Reach'
linewidth = 3
time = '60Min'

df = pd.read_csv(f, skiprows=[1])
df['Posted'] = pd.to_datetime(df['Posted'])
df = df[(df['Posted'] > '2017-08-02 12:00:00') & (df['Posted'] < '2017-08-17 12:00:00')]
df = df.set_index(['Posted'])

g30 = df.groupby(pd.TimeGrouper(freq=time))
xy = g30.aggregate(sum)[reach].cumsum()
plt.plot(xy[~np.isnan(xy)], '-', c='m', lw=linewidth, markeredgecolor='none',
         label='FB Total Post Reach')

#Twitter
dt = pd.read_csv(t)
dt['time'] = pd.to_datetime(dt['time'])
dt = dt[(dt['time'] > '2017-08-02 12:00:00') & (dt['time'] < '2017-08-17 12:00:00')]
dt = dt.set_index(['time'])

g30 = dt.groupby(pd.TimeGrouper(freq=time))
xy = g30.aggregate(sum)['impressions'].cumsum()

plt.plot(xy[~np.isnan(xy)], '-', lw=linewidth, c='c', markeredgecolor='none',
         label='Twitter Impressions')

plt.legend(loc='best')
plt.show()
